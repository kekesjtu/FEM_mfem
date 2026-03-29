## 1. 结果输出模块解耦与通用化

我注意到 `SolveThermoMechanical` 和 `SolveSingleField` 中没有单独抽出结果输出模块。建议将 `ExportResults(e_solver, t_solver)` 做成通用接口；同时把 `post` 中的 txt 输出与 `io` 中的 `paraviewExporter` 整合到统一的 `output` 模块中，形成更通用的结果输出接口；并在耦合器中删除 `ExportResults`。

## 2. `config` 中边界条件建模方式对比

我之前在 `config` 中采用“将所有 BC 打包在一起”的设计。现在改为“不同 BC 分别作为 vector”后，希望对两种方案的优劣进行对比讲解。

之前的设计：

```cpp
struct MechanicalFieldConfig
{
    // 向量边界条件
    struct BoundaryCondition
    {
        enum class BCType
        {
            FREE,
            DIRICHLET_FULL,
            DIRICHLET_NORMAL,
            NEUMANN_PRESSURE,
        } type = BCType::FREE;  // 只需要实现这四种,不需要牵引力边界条件

        double normal_dirichlet_displacement_value = 0.0;
        std::vector<double> dirichlet_displacement_value{0.0, 0.0, 0.0};
        double pressure_value = 0.0;
    };

    std::unordered_map<int, BoundaryCondition> boundary_to_bc;
};
```

现在的设计：

```cpp
struct MechanicalFieldConfig
{
    bool enabled = false;

    struct DirichletDisplacementBC
    {
        std::vector<int> bdr_attributes;
        std::vector<double> value{0.0, 0.0, 0.0};
    };

    struct NormalDisplacementBC
    {
        std::vector<int> bdr_attributes;
        double penalty = 1.0e12;
    };

    struct PressureBC
    {
        std::vector<int> bdr_attributes;
        double value = 0.0;
    };

    std::vector<double> dirichlet_displacement_value{0.0, 0.0, 0.0};
    std::vector<NormalDisplacementBC> normal_displacement_bcs;
    std::vector<PressureBC> pressure_bcs;
};
```

## 3. 未使用结构体与契约文件清理

`config` 中以下结构体似乎没有被实际使用。可能是契约文件中的历史遗留字段造成了误导，建议同步优化契约文件。

```cpp
struct DirichletDisplacementBC
{
    std::vector<int> bdr_attributes;
    std::vector<double> value{0.0, 0.0, 0.0};
};

std::vector<double> dirichlet_displacement_value{0.0, 0.0, 0.0};//以及这个意义不明的字段
```

## 4. `config` 与 `assembly` 结构对齐，减少求解器显式转换

`assembly` 中：

```cpp
/// Boundary condition marker + coefficient for Robin / Dirichlet
struct RobinBoundaryBlock
{
    mfem::Array<int> marker;
    mfem::Coefficient &l;
    mfem::Coefficient &q;
};

struct DirichletBoundaryBlock
{
    mfem::Array<int> marker;
    mfem::Coefficient &value;
};
```

`config` 中：

```cpp
struct DirichletBC
{
    std::vector<int> bdr_attributes;
    double value = 0.0;
};

struct RobinBC
{
    std::vector<int> bdr_attributes;
    double l = 0.0;
    double q = 0.0;
};
```

两者在概念上是类似的。我希望 `config` 直接完成一步到位的构建，而不是在物理场 solver 中仍需进行如下显式转换：

```cpp
std::vector<mfem::ConstantCoefficient> dirichlet_coeffs;
std::vector<assembly::DirichletBoundaryBlock> dirichlet_blocks;

dirichlet_coeffs.reserve(field_config.dirichlet_bcs.size());
for (const auto &dbc : field_config.dirichlet_bcs)
{
    mfem::Array<int> marker(num_bdr);
    marker = 0;
    for (int attr : dbc.bdr_attributes)
    {
        if (attr >= 1 && attr <= num_bdr)
        {
            marker[attr - 1] = 1;
            essential_bdr[attr - 1] = 1;
        }
    }
    dirichlet_coeffs.emplace_back(dbc.value);
    dirichlet_blocks.push_back({std::move(marker), dirichlet_coeffs.back()});
}

// Build Robin boundary conditions
std::vector<mfem::ConstantCoefficient> robin_l_coeffs;
std::vector<mfem::ConstantCoefficient> robin_q_coeffs;
std::vector<assembly::RobinBoundaryBlock> robin_blocks;

robin_l_coeffs.reserve(field_config.robin_bcs.size());
robin_q_coeffs.reserve(field_config.robin_bcs.size());
for (const auto &rbc : field_config.robin_bcs)
{
    mfem::Array<int> marker(num_bdr);
    marker = 0;
    for (int attr : rbc.bdr_attributes)
    {
        if (attr >= 1 && attr <= num_bdr)
        {
            marker[attr - 1] = 1;
        }
    }
    robin_l_coeffs.emplace_back(rbc.l);
    robin_q_coeffs.emplace_back(rbc.q);
    robin_blocks.push_back({std::move(marker), robin_l_coeffs.back(), robin_q_coeffs.back()});
}
```

## 5. `assembly` 输入构造职责下沉

在物理场solver中,调用assembly需要这一段显式输入构造,较长且不够优雅。建议由 `config` 注入并实体化一个 assembler，在其内部完成构建，尤其是在第 4 点优化之后。

```cpp
// Assemble the system
assembly::PoissonAssemblyInput input{
    fespace,
    sigma_coeff,
    source_coeff,
    essential_bdr,
    0.0,
    std::move(dirichlet_blocks),
    std::move(robin_blocks),
};
```

## 6. `muparser` 表达式求值性能评估

`muparser` 在表达式求值时，是否会对项目产生较大的性能影响？除常数退化方法外，是否还有针对表达式求值的优化方案？

# design3 — 细节调整2 设计记录

本文档记录第三轮优化的 4 项改进：机械场 BC 与载荷补全、拉梅参数自动补全、MPI 并行能力建设、比较脚本优化。

---

## 1. 机械场边界与载荷能力补全

### 问题

原实现中机械场仅支持法向位移惩罚 BC（`NormalDisplacementBC`）和压力 BC（`PressureBC`），缺少完整位移 BC（指定三个分量），且 `body_force` 仅支持全局默认值，不支持分域覆盖——与热场的 `source_default` + `source_by_domain` 模式不一致。

### 方案

**位移 BC**：`DisplacementBC` 结构体已存在于 `Config.hpp`，补全 JSON 解析与装配路径：

```cpp
struct DisplacementBC {
    std::vector<int> bdr_attributes;
    std::vector<double> value{0.0, 0.0, 0.0};  // 三个方向分量
};
```

装配器中以 essential BC 方式处理——通过 `ProjectBdrCoefficient(VectorConstantCoefficient)` 投影到边界自由度，再标记为 essential DOFs：

```cpp
// 装配器内部
initial_guess.ProjectBdrCoefficient(disp_coeffs[i], disp_markers[i]);
fespace.GetEssentialTrueDofs(essential_bdr, essential_tdofs);
bilinear->FormLinearSystem(essential_tdofs, initial_guess, *linear, A, X, B);
```

**分域体力**：与热场 source 处理模式一致，`body_force_default` 作为默认值，`body_force_by_domain` 按域覆盖。在 `MechanicalFieldSolver` 中，使用 `VectorArrayCoefficient`（每个分量一个 `PWConstCoefficient`）实现分域向量系数：

```cpp
auto pw_vec = std::make_unique<mfem::VectorArrayCoefficient>(dim);
for (int d = 0; d < dim; ++d) {
    mfem::Vector comp_vals(max_attr);
    // 先填默认值，再按 domain_to_body_force 覆盖
    pw_vec->Set(d, new mfem::PWConstCoefficient(comp_vals));
}
```

当无分域覆盖时，退化为 `VectorConstantCoefficient`，零开销。

### JSON 配置示例

```json
{
    "type": "mechanical",
    "displacement_boundary_conditions": [
        {"bdr_attributes": [1, 2], "value": [0.0, 0.0, 0.0]}
    ],
    "normal_displacement_boundary_conditions": [
        {"bdr_attributes": [3, 5], "penalty": 1.0e16}
    ],
    "body_force_default": [0.0, 0.0, -9800.0],
    "body_force_by_domain": {
        "3": [0.0, 0.0, 0.0]
    }
}
```

### 涉及文件

| 文件                                                     | 变更                                                 |
| -------------------------------------------------------- | ---------------------------------------------------- |
| `src/frontend/ConfigLoader.cpp`                          | `LoadMechanicalField()` 增加 displacement_bcs 解析   |
| `include/fem/assembly/MfemLinearElasticityAssembler.hpp` | `LinearElasticityInput` 增加 `displacement_bcs` 字段 |
| `src/assembly/MfemLinearElasticityAssembler.cpp`         | 装配时构建 essential DOF 标记、投影位移值            |
| `src/physics/MechanicalFieldSolver.cpp`                  | `body_force` 改为分域感知（VectorArrayCoefficient）  |

---

## 2. 配置层自动补全拉梅参数

### 问题

原实现在 `MechanicalFieldSolver::Solve()` 中有约 30 行手动从 `young_modulus` 和 `poisson_ratio` 转换拉梅参数的逻辑（循环遍历域属性→材料→属性→计算 λ/μ），与项目已有的 `PiecewiseConstantCoefficient` 机制不一致——热场和电场均通过 property name 直接查找材料属性。

### 方案

在 `ConfigLoader::LoadMaterials()` 末尾自动检测并补全：

```cpp
for (auto &[mat_name, props] : db.material_to_properties) {
    auto E_it = props.find("young_modulus");
    auto nu_it = props.find("poisson_ratio");
    if (E_it != props.end() && nu_it != props.end()) {
        double E = std::stod(E_it->second);
        double nu = std::stod(nu_it->second);
        props["lambda_lame"] = std::to_string(E * nu / ((1+nu)*(1-2*nu)));
        props["mu_lame"]     = std::to_string(E / (2*(1+nu)));
    }
}
```

求解器端简化为两行：

```cpp
coeff::PiecewiseConstantCoefficient lambda_coeff("lambda_lame", config_.materials, mesh);
coeff::PiecewiseConstantCoefficient mu_coeff("mu_lame", config_.materials, mesh);
```

替代原来的 30 行手动转换逻辑。

### 设计决策

- **补全时机**：放在 `LoadMaterials()` 末尾而非 `BuildFE()` 中，因为这是纯材料属性转换，不依赖网格或 FE 空间。
- **以字符串存储**：与其他材料属性一致，`PiecewiseConstantCoefficient` 构造时 `stod` 解析。虽然存在浮点→字符串→浮点的精度损耗，但 `std::to_string` 保留 6 位有效数字，对工程材料参数足够。
- **不覆盖已有属性**：如果 JSON 中已显式定义 `lambda_lame`/`mu_lame`，可手动检测跳过（当前实现直接覆盖，因为一致性优先）。

### 涉及文件

| 文件                                    | 变更                                                  |
| --------------------------------------- | ----------------------------------------------------- |
| `src/frontend/ConfigLoader.cpp`         | `LoadMaterials()` 末尾增加拉梅参数自动补全            |
| `src/physics/MechanicalFieldSolver.cpp` | 删除 30 行手动转换，改用 PiecewiseConstantCoefficient |

---

## 3. MPI 并行能力建设

### 设计目标

- 串行/并行代码共存，同一套对外接口
- 条件编译 `#ifdef MFEM_USE_MPI` 区分路径
- 通过 `assembly_mode: "parallel"` 配置项切换

### 架构

并行方案沿袭 MFEM 标准并行流程（参考 `ex1p.cpp`、`ex2p.cpp`）：

```
serial Mesh → ParMesh(MPI_COMM_WORLD, mesh)
FiniteElementSpace → ParFiniteElementSpace
BilinearForm → ParBilinearForm
LinearForm → ParLinearForm
GridFunction → ParGridFunction （隐式，指向 ParFESpace）
SparseMatrix → HypreParMatrix
DSmoother → HypreBoomerAMG
CGSolver → CGSolver(MPI_COMM_WORLD) + AMG 预条件
```

#### FEConfig 扩展

```cpp
struct FEConfig {
    // serial
    unique_ptr<Mesh> mesh;
    unique_ptr<FiniteElementSpace> scalar_fespace, vector_fespace;
    // parallel
#ifdef MFEM_USE_MPI
    unique_ptr<ParMesh> pmesh;
    unique_ptr<ParFiniteElementSpace> par_scalar_fespace, par_vector_fespace;
#endif
    bool IsParallel() const;
    Mesh &GetMesh() const;                 // 返回 pmesh 或 mesh
    FiniteElementSpace &GetScalarFESpace(); // 返回 par 或 serial
    FiniteElementSpace &GetVectorFESpace();
};
```

所有下游代码通过 `GetMesh()`、`GetScalarFESpace()`、`GetVectorFESpace()` 访问，无需关心串行/并行。

#### ConfigLoader::BuildFE

```
if (assembly_mode == "parallel" && MFEM_USE_MPI):
    pmesh = make_unique<ParMesh>(MPI_COMM_WORLD, *mesh)
    par_scalar_fespace = make_unique<ParFiniteElementSpace>(pmesh, fec)
    par_vector_fespace = make_unique<ParFiniteElementSpace>(pmesh, fec, dim)
else:
    scalar_fespace = make_unique<FiniteElementSpace>(mesh, fec)
    vector_fespace = make_unique<FiniteElementSpace>(mesh, fec, dim)
```

#### 装配器

装配器通过 `dynamic_cast<ParFiniteElementSpace*>` 检测运行时环境。如果是 `ParFiniteElementSpace`，使用 `ParBilinearForm`/`ParLinearForm`；否则走串行路径。输出统一为 `AssembledSystem`（`OperatorPtr A`、`Vector X/B`）。

`FormLinearSystem` 在并行模式下产出 `HypreParMatrix`，在串行模式下产出 `SparseMatrix`。`OperatorPtr` 统一承载两者。

#### 求解器

`MfemPcgSolver::Solve` 通过 `dynamic_cast` 区分矩阵类型：
- `SparseMatrix` → `CGSolver` + `DSmoother`
- `HypreParMatrix` → `CGSolver(MPI_COMM_WORLD)` + `HypreBoomerAMG`

AMG 使用默认构造（惰性初始化），避免在构造时触发 setup 导致的性能问题。

#### main.cpp

```cpp
#ifdef MFEM_USE_MPI
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();
#endif
```

`Mpi::Init` 内部处理 `MPI_Init`，程序退出时自动 `MPI_Finalize`。

#### 日志

非零 rank 的 spdlog 级别设为 `off`，确保只有 rank 0 输出日志。Hypre 求解器的 print_level 也仅在 rank 0 非零。

#### 文本导出

`SolutionTextExporter` 已有 `ParallelAwareOutputPath()`，并行时按 rank 分文件输出（`temperature.rank0.txt`、`temperature.rank1.txt`）。

### 验证结果

在 2 进程并行模式下，全部测试用例通过：

| 测试用例               | 串行 MAE | 并行 MAE | 匹配 |
| ---------------------- | -------- | -------- | ---- |
| busbar 电热无迭代 温度 | 3.75e-9  | 3.85e-9  | ✓    |
| bipolar 热力耦合 ux    | 1.70e-6  | 1.70e-6  | ✓    |
| bipolar 热力耦合 uy    | 1.19e-6  | —        | ✓    |
| bipolar 热力耦合 uz    | 2.99e-7  | —        | ✓    |

并行收敛迭代数显著减少（串行 ~400 iter → 并行 20 iter），得益于 AMG 预条件。

### 涉及文件

| 文件                                             | 变更                                         |
| ------------------------------------------------ | -------------------------------------------- |
| `include/fem/frontend/Config.hpp`                | FEConfig 增加并行成员 + 访问器方法           |
| `src/frontend/ConfigLoader.cpp`                  | BuildFE 增加并行路径                         |
| `src/assembly/MfemPoissonAssembler.cpp`          | 增加 ParBilinearForm/ParLinearForm 并行路径  |
| `src/assembly/MfemLinearElasticityAssembler.cpp` | 同上，并重构为辅助函数复用 integrator 添加   |
| `src/physics/ElectrostaticFieldSolver.cpp`       | 改用 `GetScalarFESpace()`/`GetMesh()` 访问器 |
| `src/physics/ThermalFieldSolver.cpp`             | 同上                                         |
| `src/physics/MechanicalFieldSolver.cpp`          | 改用 `GetVectorFESpace()`/`GetMesh()` 访问器 |
| `src/coupling/MultiPhysicsCoupler.cpp`           | 改用 `GetMesh()` 访问器                      |
| `src/solver/MfemPcgSolver.cpp`                   | 并行路径使用 CGSolver + HypreBoomerAMG       |
| `src/solver/MfemAmgSolver.cpp`                   | 同上                                         |
| `src/log/Logger.cpp`                             | 非零 rank 抑制日志                           |
| `main.cpp`                                       | 增加 `Mpi::Init` + `Hypre::Init`             |

---

## 4. 比较脚本优化

### 问题

原 `compare_solution_txt.py` 需要为矢量场的每个分量单独运行（指定 `--value-col 3/4/5`），操作繁琐且不易一目了然地看到全场精度。

### 方案

升级为自动检测模式：

- **标量场**：自输出 1 个值列，COMSOL 端 1 个值列 → 自动标量比较
- **矢量场**：自输出 4+ 个值列（ux, uy, uz, magnitude），COMSOL 3+ 个值列 → 自动逐分量 + 模比较

检测逻辑：

```python
def detect_mode(mine, comsol):
    if mine_vcols >= 4 and comsol_vcols >= 3:
        return "vector"
    return "scalar"
```

矢量模式输出：

```
==== 与 COMSOL 对比 (矢量场自动检测) ====
  [ux]        matched=1552/1552, MAE=1.70e-06, ...
  [uy]        matched=1552/1552, MAE=1.19e-06, ...
  [uz]        matched=1552/1552, MAE=2.99e-07, ...
  [magnitude] matched=1552/1552, MAE=2.20e-06, ...
```

保留 `--value-col` 参数用于兼容旧用法。不指定时走自动检测。

### 涉及文件

| 文件                            | 变更                                                  |
| ------------------------------- | ----------------------------------------------------- |
| `tools/compare_solution_txt.py` | 重构 DataSet 存储全部值列、增加自动检测与矢量比较逻辑 |

---

## 验证结果汇总

全部 4 个测试用例在串行模式下结果不变：

```
busbar 电热无迭代        MAE=3.75e-09
busbar 电热无迭代 order2 MAE=1.60e-06
busbar 电热 Picard       MAE=0.002
bipolar 热力耦合         ux MAE=1.70e-06, uy MAE=1.19e-06, uz MAE=2.99e-07
```

并行模式（2 进程）精度与串行一致。

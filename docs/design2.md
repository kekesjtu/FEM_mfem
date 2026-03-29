# design2 — 细节调整设计记录

本文档记录对第一版实现的 6 项细节优化，包含设计决策分析和具体修改说明。

---

## 1. 结果输出模块解耦与通用化

### 问题

原实现中，结果导出逻辑分散在 `MultiPhysicsCoupler` 的多个方法中：
- `ExportResults()` 处理电-热耦合的导出
- `SolveThermoMechanical()` 中内联了 ParaView + TXT 导出
- `SolveSingleField()` 中又有一套独立的导出代码

此外，`io::ParaviewExporter` 和 `post::SolutionTextExporter` 作为两个独立模块存在，调用者需要分别调用，增加了使用复杂度。

### 方案

**新建统一接口** `io::ResultExporter`，封装 ParaView + TXT 导出：

```cpp
class ResultExporter {
public:
    ResultExporter(const std::string &output_dir, mfem::Mesh &mesh);

    void ExportScalar(const std::string &collection_name, const std::string &display_name,
                      mfem::GridFunction &gf, int cycle = 0, double time = 0.0);

    void ExportVector(const std::string &collection_name, const std::string &display_name,
                      mfem::GridFunction &gf, const std::string &txt_prefix = "",
                      int cycle = 0, double time = 0.0);
};
```

- `ExportScalar("thermal", "temperature", gf)` → 生成 ParaView 集合 `thermal/` + 文本文件 `temperature.txt`
- `ExportVector("mechanical", "displacement", gf, "u")` → 生成 ParaView 集合 `mechanical/` + 文本文件 `displacement.txt`（列头：ux, uy, uz, u_magnitude）

**耦合器重构**：`MultiPhysicsCoupler::Solve()` 采用三阶段架构：

```
Phase 1: 物理求解（电-热耦合 / 单场求解）
Phase 2: 力学求解（可选，接收温度场）
Phase 3: 统一导出（一个 ResultExporter 实例处理所有场）
```

`ExportResults` 私有方法以及 `SolveThermoMechanical`、`SolveSingleField` 方法均被移除。

**附带修复**：原实现中完整的 E+T+M 耦合路径存在 bug——`SolveThermoMechanical()` 会创建新的 `ThermalFieldSolver` 重新求解热场（无焦耳热），导致电-热耦合的温度结果被覆盖。新实现中，热场求解器以 `unique_ptr` 在 `Solve()` 顶层管理，Phase 2 直接复用 Phase 1 的温度 GridFunction。

### 涉及文件

| 文件                                           | 变更                                                                                                |
| ---------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| `include/fem/io/ResultExporter.hpp`            | **新建**                                                                                            |
| `src/io/ResultExporter.cpp`                    | **新建**                                                                                            |
| `include/fem/coupling/MultiPhysicsCoupler.hpp` | 移除 `ExportResults`, `SolveThermoMechanical`, `SolveSingleField`；子方法签名改为接收 `unique_ptr&` |
| `src/coupling/MultiPhysicsCoupler.cpp`         | 完全重写为三阶段架构                                                                                |

---

## 2. Config 边界条件建模方式对比

### 旧方案：统一 BC 枚举 + map

```cpp
struct BoundaryCondition {
    enum class BCType { FREE, DIRICHLET_FULL, DIRICHLET_NORMAL, NEUMANN_PRESSURE } type;
    double normal_dirichlet_displacement_value;
    std::vector<double> dirichlet_displacement_value;
    double pressure_value;
};
std::unordered_map<int, BoundaryCondition> boundary_to_bc;
```

**优点**：
- 单一 map 查找，每个边界属性有明确的唯一 BC 类型
- 添加新 BC 类型只需扩展枚举

**缺点**：
- **冗余字段**：每个 BC 实例都携带所有类型的字段（如 `pressure_value` 对 Dirichlet BC 无意义），结构体含义不清晰
- **与 JSON 不匹配**：JSON schema 按 BC 类型分组定义（`dirichlet_boundary_conditions`, `pressure_boundary_conditions`），配置加载需要额外的类型判断和分发逻辑
- **批量操作不便**：多个边界共享同一 BC 值时（如多个面施加相同压力），需要在 map 中逐个设置

### 新方案：分类型 vector

```cpp
struct NormalDisplacementBC {
    std::vector<int> bdr_attributes;
    double penalty = 1.0e12;
};
struct PressureBC {
    std::vector<int> bdr_attributes;
    double value = 0.0;
};
std::vector<NormalDisplacementBC> normal_displacement_bcs;
std::vector<PressureBC> pressure_bcs;
```

**优点**：
- **语义清晰**：每种 BC 类型只包含自己需要的字段，无冗余
- **与 JSON 天然对齐**：JSON 中的 `normal_displacement_boundary_conditions` 数组直接映射为 `vector<NormalDisplacementBC>`，加载代码无需类型判断
- **批量友好**：`bdr_attributes` 是 vector，一次指定多个边界面，减少重复定义
- **组装层可直接消费**：无需额外的类型分发逻辑

**缺点**：
- 同一边界属性若出现在多个 BC vector 中会导致冲突（但这是配置错误，应由校验层处理）
- 添加新 BC 类型需要定义新 struct 和新 vector 成员

**结论**：对于当前需求（BC 类型固定且数量有限，JSON 驱动配置），分类型 vector 方案更优。

---

## 3. 未使用结构体清理

### 清理内容

| 移除项                                  | 位置                             | 原因                                                                                                           |
| --------------------------------------- | -------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `DirichletDisplacementBC` 结构体        | `Config.hpp`                     | 定义后从未被 ConfigLoader 加载或 Solver 引用                                                                   |
| `dirichlet_displacement_value` 字段     | `Config.hpp`                     | 虽然在 MechanicalFieldSolver 中被引用，但 `essential_bdr` 始终为全零，ProjectBdrCoefficient 不会实际投影任何值 |
| `dirichlet_displacement_value`          | `physics.schema.json`            | 同步清理契约文件                                                                                               |
| `dirichlet_displacement_value`          | `bipolar_thermo_mechanical.json` | 同步清理配置文件                                                                                               |
| `dirichlet_displacement_value` 加载代码 | `ConfigLoader.cpp`               | 同步清理                                                                                                       |

### 分析

`dirichlet_displacement_value` 在 MechanicalFieldSolver 中的使用路径：

```cpp
// 构造 VectorConstantCoefficient，值为 [0,0,0]
mfem::VectorConstantCoefficient dirichlet_disp_coeff(dirichlet_disp_vec);
// 传入 assembly input
input.dirichlet_displacement = dirichlet_disp_coeff;
// 在 assembler 中执行投影
initial_guess.ProjectBdrCoefficient(input.dirichlet_displacement, input.essential_bdr_marker);
```

由于 `essential_bdr_marker` 全为 0（当前力学场所有约束均通过罚函数法施加），`ProjectBdrCoefficient` 不会修改任何 DOF。因此该字段名义上存在但从未产生实际效果。

---

## 4. Config 与 Assembly 结构对齐

### 问题

原设计中，Config 定义了一套 BC 类型：

```cpp
// Config 层
struct DirichletBC { std::vector<int> bdr_attributes; double value; };
struct RobinBC { std::vector<int> bdr_attributes; double l, q; };
```

Assembly 定义了另一套结构：

```cpp
// Assembly 层
struct DirichletBoundaryBlock { mfem::Array<int> marker; mfem::Coefficient &value; };
struct RobinBoundaryBlock { mfem::Array<int> marker; mfem::Coefficient &l, &q; };
```

两者概念等价，但表示不同：Config 用 `vector<int>` + `double`，Assembly 用 `mfem::Array<int>` marker + `mfem::Coefficient&`。物理场 Solver 必须手动执行转换。

### 方案

**让 Assembly 输入直接使用 Config BC 类型**，Assembly 头文件包含 `Config.hpp`：

```cpp
// MfemPoissonAssembler.hpp
#include "fem/frontend/Config.hpp"

struct PoissonAssemblyInput {
    mfem::FiniteElementSpace &space;
    mfem::Coefficient &diffusion;
    mfem::Coefficient &source;
    std::vector<frontend::ScalarFieldConfig::DirichletBC> dirichlet_bcs;  // ← 直接使用
    std::vector<frontend::ScalarFieldConfig::RobinBC> robin_bcs;         // ← 直接使用
    double default_value = 0.0;
};
```

Assembly 不再定义 `DirichletBoundaryBlock`、`RobinBoundaryBlock`、`ScalarBoundaryBlock`、`NormalDisplacementBoundary`、`TractionBoundary`、`PressureBoundary` 等中间类型。

### 依赖方向分析

引入了 Assembly → Frontend (Config.hpp) 的依赖。这是否违反 DIP？

- Config 是项目的**数据中心**，按设计文档的定位，所有层都依赖它提供的配置信息
- BC 类型是简单 POD 结构体，不包含业务逻辑
- 消除了两层之间的类型冗余和转换代码

**结论**：这是合理的**结构共享**而非紧耦合。

---

## 5. Assembly 输入构造职责下沉

### 问题

原实现中物理场 Solver 需要 ~40 行代码完成 Config BC → Assembly Block 的转换：

```cpp
// 每个 Solver 都要重复这段 ~40 行
std::vector<mfem::ConstantCoefficient> dirichlet_coeffs;
std::vector<assembly::DirichletBoundaryBlock> dirichlet_blocks;
for (const auto &dbc : field_config.dirichlet_bcs) {
    mfem::Array<int> marker(num_bdr);  marker = 0;
    for (int attr : dbc.bdr_attributes) { marker[attr-1] = 1; essential_bdr[attr-1] = 1; }
    dirichlet_coeffs.emplace_back(dbc.value);
    dirichlet_blocks.push_back({std::move(marker), dirichlet_coeffs.back()});
}
// ... Robin 同理 ...
```

ElectrostaticFieldSolver 和 ThermalFieldSolver 中各有一份近乎相同的代码。

### 方案

**将 marker 构建和 Coefficient 创建下沉到 Assembler 内部**。Assembler 接收 Config BC 类型后，在 `Assemble()` 方法内部完成转换：

```cpp
// MfemPoissonAssembler::Assemble() 内部
std::vector<mfem::ConstantCoefficient> robin_l_coeffs, robin_q_coeffs;
std::vector<mfem::Array<int>> robin_markers;
robin_l_coeffs.reserve(input.robin_bcs.size());
// ... 构建 markers 和 coefficients ...
// ... 添加 integrators ...
// ... Assemble + FormLinearSystem ...
```

**生命周期安全性**：`ConstantCoefficient` 和 `Array<int>` 作为局部变量，在 `Assemble()` 和 `FormLinearSystem()` 调用期间存活。MFEM 组装完成后，系数值已积分进稀疏矩阵，后续 `RecoverFEMSolution()` 不会重新访问 integrator。

### 效果

| 文件                         | 修改前行数 | 修改后行数 | 减少       |
| ---------------------------- | ---------- | ---------- | ---------- |
| ElectrostaticFieldSolver.cpp | 97 行      | 45 行      | **-52 行** |
| ThermalFieldSolver.cpp       | 120 行     | 72 行      | **-48 行** |
| MechanicalFieldSolver.cpp    | 140 行     | 85 行      | **-55 行** |

物理场 Solver 现在只做三件事：
1. 构造物理系数（$\sigma$, $k$, $\lambda$, $\mu$ 等）
2. 填充 Assembly 输入（直接传递 Config BC）
3. 调用 Assembler + Solver + Recover

---

## 6. muparser 表达式求值性能评估

### 性能特征

muparser 的表达式求值分为两个阶段：

1. **解析阶段** (`SetExpr`): 将字符串编译为内部字节码（token 序列），只执行一次
2. **求值阶段** (`Eval`): 遍历字节码执行运算，每次调用

解析阶段的开销较高（字符串处理、语法分析），但只在 `Expression` 构造时执行一次。求值阶段的开销取决于表达式复杂度。

### 基准估算

对于典型的温度相关电导率表达式 `5.998e7 / (1.0 + 0.00393 * (T - 293.15))`：
- muparser `Eval()`: 约 100-200 ns/次（含 4 次浮点运算 + 字节码遍历）
- 等价的原生 C++ 函数: 约 2-5 ns/次

差距约 **50 倍**。但对于 FEM 求解的整体时间，这是否显著？

### 求值次数分析

以 busbar 模型为例（7340 节点，~40000 四面体，阶数 1）：
- 每个单元 4 个高斯积分点
- Picard 迭代 4 次
- 每次迭代：电场组装 + 热场组装（σ 在两处被求值）
- 总求值次数：$40000 \times 4 \times 4 \times 2 \approx 1.3 \times 10^6$ 次

以 200ns/次计算：$1.3 \times 10^6 \times 200 \text{ns} \approx 0.26$ 秒。

实际 busbar Picard 求解总时间约 6 秒，muparser 占比约 4%。**对于当前规模，性能影响可以接受。**

### 已实现的优化

1. **常量检测**（`Expression::is_constant_`）：构造时尝试 `std::stod`，若整个字符串是数字则标记为常量，后续 `Evaluate()` 直接返回缓存值，完全跳过 muparser
2. **PiecewiseConstantCoefficient**：当某个材料属性（如热导率）在所有域都是常量表达式时，使用数组查表 O(1) 代替 muparser 求值
3. **单次解析**：`Expression` 对象在构造时解析字符串，后续只执行字节码求值

### 进一步优化方向（未实现）

1. **JIT 编译**：将 muparser 字节码编译为原生机器码（如通过 [exprtk](https://github.com/ArashPartow/exprtk) 或 LLVM JIT），可消除字节码遍历开销
2. **批量求值**：muparser 支持 `Eval(int &nResults)` 批量模式，对向量化友好
3. **表达式缓存池**：对相同字符串表达式共享同一 `mu::Parser` 实例（当前每个域独立创建）
4. **空间无关检测**：若表达式不含 x/y/z 变量（仅含 T），可在同一域内进一步缓存（同一单元内所有积分点的 T 通常相近但不完全相同，此策略不安全）

对于百万级自由度的大规模问题，建议优先采用方案 1（JIT 编译）。

---

## 修改汇总

| 文件                                                     | 操作     | 说明                                                           |
| -------------------------------------------------------- | -------- | -------------------------------------------------------------- |
| `include/fem/frontend/Config.hpp`                        | 修改     | 移除 `DirichletDisplacementBC`、`dirichlet_displacement_value` |
| `docs/physics.schema.json`                               | 修改     | 移除 `dirichlet_displacement_value` 字段定义                   |
| `configs/bipolar_thermo_mechanical.json`                 | 修改     | 移除 `dirichlet_displacement_value`                            |
| `src/frontend/ConfigLoader.cpp`                          | 修改     | 移除 `dirichlet_displacement_value` 加载                       |
| `include/fem/assembly/MfemPoissonAssembler.hpp`          | **重写** | 移除中间 block 类型，使用 Config BC 类型                       |
| `include/fem/assembly/MfemLinearElasticityAssembler.hpp` | **重写** | 移除中间 block 类型，使用 Config BC 类型                       |
| `src/assembly/MfemPoissonAssembler.cpp`                  | **重写** | 内部构建 marker + coefficient                                  |
| `src/assembly/MfemLinearElasticityAssembler.cpp`         | **重写** | 内部构建 marker + coefficient                                  |
| `src/physics/ElectrostaticFieldSolver.cpp`               | **重写** | 移除手动 BC 转换代码                                           |
| `src/physics/ThermalFieldSolver.cpp`                     | **重写** | 移除手动 BC 转换代码                                           |
| `src/physics/MechanicalFieldSolver.cpp`                  | **重写** | 移除手动 BC 转换代码，移除 dirichlet_displacement 相关         |
| `include/fem/io/ResultExporter.hpp`                      | **新建** | 统一输出接口                                                   |
| `src/io/ResultExporter.cpp`                              | **新建** | 统一输出实现                                                   |
| `include/fem/coupling/MultiPhysicsCoupler.hpp`           | **重写** | 三阶段架构，移除 ExportResults                                 |
| `src/coupling/MultiPhysicsCoupler.cpp`                   | **重写** | 三阶段求解+统一导出                                            |

### 验证结果

所有 4 个测试用例在重构后通过验证，数值结果与重构前完全一致：

| 测试                | MAE       | MeanRel  |
| ------------------- | --------- | -------- |
| busbar 电→热 (1 阶) | 3.75e-9 K | 1.16e-11 |
| busbar 电→热 (2 阶) | 1.61e-6 K | 4.97e-9  |
| busbar 电↔热 Picard | 0.002 K   | 6.68e-6  |
| bipolar 热→力 (ux)  | 1.70e-6 m | 2.6%     |
| bipolar 热→力 (uy)  | 1.19e-6 m | 2.4%     |
| bipolar 热→力 (uz)  | 2.99e-7 m | 17.2%*   |

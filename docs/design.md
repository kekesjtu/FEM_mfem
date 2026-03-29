# 多物理场有限元求解器框架 — 设计文档

## 1. 概览

本项目实现了一个基于 **MFEM** 的多物理场有限元求解器框架，支持：

- **静电场** (Electrostatic)：求解电势分布 $V$
- **热场** (Thermal)：求解温度分布 $T$，可接收焦耳热源
- **力学场** (Mechanical)：求解位移分布 $\mathbf{u}$，可接收热膨胀载荷

以及以下耦合策略：

| 耦合模式                | 描述                                    |
| ----------------------- | --------------------------------------- |
| 电→热 单向耦合          | 电场求解 → 焦耳热 $Q=\sigma             | \nabla V | ^2$ → 热场求解 |
| 电↔热 双向耦合 (Picard) | 电导率 $\sigma(T)$ 依赖温度，迭代至收敛 |
| 热→力 单向耦合          | 温度场 → 热膨胀应变 → 力学求解          |
| 电↔热 + 热→力           | 完整热-力-电多场耦合                    |

---

## 2. 系统架构

```
┌─────────────────────────────────────────────────────────┐
│                    main.cpp                              │
│              Application::Run(configPath)                │
└──────────────────────┬──────────────────────────────────┘
                       │
          ┌────────────▼────────────┐
          │   MultiPhysicsCoupler   │  ← 耦合层 (coupling/)
          │  协调求解顺序与数据交换   │
          └────┬───────┬───────┬────┘
               │       │       │
    ┌──────────▼──┐ ┌──▼──────┐ ┌──▼──────────────┐
    │Electrostatic│ │ Thermal │ │   Mechanical     │  ← 物理场层 (physics/)
    │FieldSolver  │ │FieldSolver│ │FieldSolver      │
    └──────┬──────┘ └──┬──────┘ └──┬──────────────┘
           │           │           │
    ┌──────▼───────────▼───────────▼──────┐
    │         Assembly Layer               │  ← 组装层 (assembly/)
    │  PoissonAssembler / ElasticityAssembler│
    └──────────────────┬──────────────────┘
                       │
    ┌──────────────────▼──────────────────┐
    │          Solver Layer                │  ← 求解层 (solver/)
    │   SolverFactory → PCG/UMFPACK/...   │
    └──────────────────┬──────────────────┘
                       │
    ┌──────────────────▼──────────────────┐
    │         Output Layer (io/)            │  ← 输出层
    │  ParaViewExporter / NodalTxtExporter │
    └─────────────────────────────────────┘
```

**数据流**：`JSON → ConfigLoader → ProjectConfig → PhysicsSolvers → Assemblers → Solvers → Export`

---

## 3. 核心设计原则

### 3.1 依赖倒置 (DIP)

- 物理场求解器**不直接依赖**组装层的具体实现。组装层定义独立的输入结构体（`PoissonAssemblyInput`、`LinearElasticityInput`），物理场求解器只需填充这些结构体。
- 组装层不依赖 `ProjectConfig`，只接收自己需要的输入数据。
- 求解器层通过工厂模式（`SolverFactory`）创建，上层不关心具体求解器类型。

### 3.2 Config 作为唯一数据中心

`ProjectConfig` 是整个系统的**单一数据源 (Single Source of Truth)**。它包含：

```cpp
struct ProjectConfig {
    SimulationConfig simulation;  // 求解器配置、输出路径
    MaterialDatabase materials;   // 材料属性（表达式或常数）
    FEConfig fe;                  // 网格、有限元空间（unique_ptr 管理）
    // 各物理场配置
    std::optional<ScalarFieldConfig> electrostatic;
    std::optional<ScalarFieldConfig> thermal;
    std::optional<MechanicalFieldConfig> mechanical;
};
```

所有下游模块从 `ProjectConfig` 获取信息，不构建额外的全局状态。

### 3.3 前后端解耦的材料系统

**前端**（JSON）定义材料属性为字符串表达式：
```json
{
  "sigma": "5.998e7 / (1.0 + 0.00393 * (T - 293.15))"
}
```

**转换链**：
1. `ConfigLoader` 将 JSON 字符串存入 `MaterialDatabase` 的 `std::unordered_map<std::string, std::string>` 属性映射
2. 物理场求解器调用 `CoefficientManager` 中的工具函数，将字符串表达式转换为 MFEM `Coefficient` 对象
3. 三种系数类型：
   - `ExpressionCoefficient`：通过 muparser 解析表达式，支持 T 变量耦合
   - `PiecewiseConstantCoefficient`：常数表达式优化（一次解析，后续直接查表）
   - `JouleHeatingCoefficient`：$Q = \sigma |\nabla V|^2$，在高斯积分点精度计算

---

## 4. 模块详细设计

### 4.1 前端层 (frontend/)

#### Expression (Expression.hpp/cpp)

封装 muparser 表达式求值引擎。

- **变量绑定**：`x, y, z, t, T`（空间坐标、时间、温度）
- **性能优化**：
  - 常量检测：`IsConstant()` 方法在初始化时尝试多组随机输入，若结果不变则标记为常量
  - 常量表达式只求值一次，后续直接返回缓存值
  - 避免对每个高斯积分点重复解析字符串

```cpp
class Expression {
    mu::Parser parser_;
    double x_, y_, z_, t_, T_;  // 绑定变量
    bool is_constant_;
    double cached_value_;
public:
    double Evaluate(double x, double y, double z, double t = 0, double T = 0);
    bool IsConstant() const;
};
```

#### ConfigLoader (ConfigLoader.hpp/cpp)

- `Load(const std::string& path) → ProjectConfig`
- 内部方法：
  - `LoadSimulation()` → 求解器类型/精度/最大迭代/输出目录
  - `LoadMaterials()` → 材料名→属性映射
  - `LoadFields()` → 电/热/力场配置（边界条件、源项、材料分配）
  - `BuildFE()` → 加载网格、创建有限元空间（标量H1/矢量H1）

#### Config.hpp 关键数据结构

**边界条件类型**：

```cpp
struct DirichletBC { std::vector<int> boundaries; double value; };
struct RobinBC { std::vector<int> boundaries; double alpha; double beta; };
struct NormalDisplacementBC { std::vector<int> boundaries; double value; double penalty; };
struct PressureBC { std::vector<int> boundaries; double value; };
```

**场配置**：

```cpp
struct ScalarFieldConfig {
    bool enabled = false;
    std::string name;
    std::unordered_map<int, std::string> domain_materials;  // domain_id → material_name
    std::unordered_map<int, double> domain_sources;         // domain_id → source_value
    std::vector<DirichletBC> dirichlet_bcs;
    std::vector<RobinBC> robin_bcs;
};

struct MechanicalFieldConfig {
    bool enabled = false;
    std::string name;
    std::unordered_map<int, std::string> domain_materials;
    std::vector<double> body_force;
    std::vector<DirichletBC> fixed_bcs;
    std::vector<NormalDisplacementBC> normal_displacement_bcs;
    std::vector<PressureBC> pressure_bcs;
    ThermalExpansionConfig thermal_expansion;
};
```

### 4.2 系数管理层 (coeff/)

#### ExpressionCoefficient

继承 `mfem::Coefficient`，在 `Eval()` 方法中：
1. 获取积分点的物理坐标 $(x, y, z)$
2. 确定所在域（通过 `T->GetAttribute()`）获取对应材料表达式
3. 若表达式依赖温度且提供了温度 GridFunction，在积分点处插值温度
4. 调用 `Expression::Evaluate()` 计算系数值

**缓存策略**：每个域-材料组合只创建一个 `Expression` 对象，避免重复解析。

#### PiecewiseConstantCoefficient

对所有表达式均为常数的场景优化：
- 初始化时一次性求值所有域的表达式
- `Eval()` 中直接按 `attribute` 查数组，零解析开销

#### JouleHeatingCoefficient

$Q = \sigma |\nabla V|^2$

- 在每个积分点：
  1. 计算 $\nabla V$（通过 `voltage_gf_->GetGradient()`）
  2. 计算当前点的 $\sigma$（通过 `sigma_coeff_->Eval()`）
  3. 返回 $\sigma \cdot |\nabla V|^2$

全部在高斯积分点精度完成，无任何节点投影。

### 4.3 组装层 (assembly/)

#### PoissonAssembler

求解标量 Poisson 类方程：$-\nabla \cdot (k \nabla u) = f$

**输入结构体**（与 Config 解耦）：

```cpp
struct PoissonAssemblyInput {
    mfem::FiniteElementSpace* fespace;
    mfem::Coefficient* diffusion_coeff;       // k
    mfem::Coefficient* source_coeff;          // f（可选）
    std::vector<ScalarBoundaryBlock> scalar_boundaries;  // Neumann
    std::vector<RobinBoundaryBlock> robin_boundaries;    // Robin
    std::vector<DirichletBoundaryBlock> dirichlet_boundaries;  // Dirichlet
};
```

**组装流程**：
1. 创建 `BilinearForm`，添加 `DiffusionIntegrator(k)`
2. 对 Robin 边界添加 `BoundaryMassIntegrator(alpha)`
3. 创建 `LinearForm`，添加 `DomainLFIntegrator(f)`
4. 对 Robin 边界添加 `BoundaryLFIntegrator(beta)`
5. 施加 Dirichlet BC（通过 `mfem::Array<int>` 标记 essential DOFs）
6. 形成 `SparseMatrix` + `Vector` 线性系统

#### LinearElasticityAssembler

求解线弹性方程：$-\nabla \cdot \boldsymbol{\sigma} = \mathbf{f}$

其中 $\boldsymbol{\sigma} = \lambda (\nabla \cdot \mathbf{u}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon}$

支持：
- 体力、压力边界、法向位移约束（罚函数法）
- 热膨胀应变：$\boldsymbol{\varepsilon}_{th} = \alpha_{th} (T - T_{ref}) \mathbf{I}$

### 4.4 求解层 (solver/)

`SolverFactory` 工厂模式：

| 求解器类型  | 实现                                        |
| ----------- | ------------------------------------------- |
| `"pcg"`     | `mfem::CGSolver` + `GSSmoother` 预条件      |
| `"umfpack"` | `mfem::UMFPackSolver`（直接法，适合小规模） |
| `"amg"`     | `mfem::CGSolver` + `HypreBoomerAMG` 预条件  |

### 4.5 物理场层 (physics/)

每个物理场求解器职责：
1. 拥有并管理自己的 `GridFunction`（解向量）
2. 从 `ProjectConfig` 读取配置
3. 构造 MFEM `Coefficient` 对象（通过 CoefficientManager）
4. 填充组装层输入结构体
5. 调用 Assembler 和 Solver 完成求解
6. 提供外部访问接口（如 `GetSolution()`）供耦合层使用

#### ElectrostaticFieldSolver

- 方程：$-\nabla \cdot (\sigma \nabla V) = 0$
- 电导率 $\sigma$ 通过 `ExpressionCoefficient` 实现，支持温度依赖
- 可接收外部 `temperature_gf` 用于 Picard 耦合

#### ThermalFieldSolver

- 方程：$-\nabla \cdot (k \nabla T) = Q$
- 热导率 $k$ 通过 `PiecewiseConstantCoefficient` 实现（通常为常数）
- 源项 $Q$ 可以是：
  - 配置文件中的分域常数源
  - `JouleHeatingCoefficient`（焦耳热耦合）
  - 两者相加
- 支持 Robin BC 实现对流换热：$q = \alpha(T - T_{ext})$

#### MechanicalFieldSolver

- 方程：$-\nabla \cdot [\lambda (\nabla \cdot \mathbf{u})\mathbf{I} + 2\mu \boldsymbol{\varepsilon}(\mathbf{u})] = \mathbf{f} + \mathbf{f}_{th}$
- Lamé 参数从 Young's modulus $E$ 和 Poisson's ratio $\nu$ 计算
- 热膨胀：温度差 $\Delta T = T - T_{ref}$ GridFunction 传入 `ThermalStrainIntegrator`
- 法向位移约束通过 `NormalDisplacementPenaltyIntegrator` 实现（罚函数法）

### 4.6 耦合层 (coupling/)

`MultiPhysicsCoupler` 根据 `ProjectConfig` 自动判断耦合策略：

```cpp
void MultiPhysicsCoupler::Solve() {
    bool has_e = config_.electrostatic.has_value() && config_.electrostatic->enabled;
    bool has_t = config_.thermal.has_value() && config_.thermal->enabled;
    bool has_m = config_.mechanical.has_value() && config_.mechanical->enabled;
    bool picard = config_.NeedsPicardIteration();

    if (has_e && has_t && picard)
        SolveElectroThermalPicard();      // 电↔热双向
    else if (has_e && has_t)
        SolveElectroThermalOneWay();      // 电→热单向
    // ... 更多组合
}
```

#### Picard 迭代 (电↔热双向耦合)

```
初始化 T = T_init (来自 Dirichlet BC 均值)

repeat:
    用当前 T 更新 σ(T)
    求解电场: -∇·(σ(T)∇V) = 0
    计算焦耳热: Q = σ(T)|∇V|²
    求解热场: -∇·(k∇T_new) = Q
    
    松弛更新: T = ω·T_new + (1-ω)·T_old  (ω = 0.7)
    
    检查收敛: ‖T_new - T_old‖/‖T_new‖ < tolerance
until converged or max_iterations reached
```

松弛因子 $\omega = 0.7$ 保证收敛稳定性。

### 4.7 输出层 (io/)

- `ParaViewExporter`：使用 MFEM 的 `ParaViewDataCollection` 导出 `.vtu` 文件
- `NodalTxtExporter`：导出节点坐标+值的文本文件，用于与 COMSOL 对比
  - 标量场：`x y z value`
  - 矢量场：`x y z ux uy uz magnitude`

---

## 5. 性能优化策略

### 5.1 Expression 常量检测

```
表达式 "5.998e7" → IsConstant() = true → 只求值一次
表达式 "5.998e7/(1+0.00393*(T-293.15))" → IsConstant() = false → 每次求值
```

对于全常量场（如热导率），自动使用 `PiecewiseConstantCoefficient`（数组查表 O(1)），避免 muparser 的字符串解析。

### 5.2 PiecewiseConstantCoefficient

预计算所有域的常数值，存入 `std::vector<double>`，`Eval()` 中只做一次数组索引：

```cpp
double Eval(ElementTransformation& T, const IntegrationPoint&) override {
    return values_[T.Attribute - 1];
}
```

### 5.3 GridFunction 梯度复用

`JouleHeatingCoefficient` 中每个积分点直接调用 `GetGradient()`，无需额外投影步骤。

---

## 6. 边界条件实现

| BC 类型   | 实现方式                                                                  | 适用场           |
| --------- | ------------------------------------------------------------------------- | ---------------- |
| Dirichlet | Essential BC，标记 DOFs 后从系统中消去                                    | 电场、热场、力学 |
| Robin     | $\alpha u + \beta$ 添加 `BoundaryMassIntegrator` + `BoundaryLFIntegrator` | 热场（对流）     |
| 法向位移  | `NormalDisplacementPenaltyIntegrator`（罚函数）                           | 力学场           |
| 压力      | `BoundaryLFIntegrator` 与法向矢量                                         | 力学场           |
| 固定约束  | Essential BC，所有分量设为 0                                              | 力学场           |

---

## 7. 文件组织

```
include/fem/
├── app/           Application.hpp          入口类
├── assembly/      MfemPoissonAssembler.hpp  标量场组装
│                  MfemLinearElasticityAssembler.hpp  力学场组装
├── coeff/         Coeffmanagaer.hpp         系数管理
├── coupling/      MultiPhysicsCoupler.hpp   耦合策略
├── frontend/      Config.hpp               数据结构定义
│                  ConfigLoader.hpp          JSON加载
│                  Expression.hpp            表达式引擎
├── integrator/    (自定义积分器)
├── io/            ParaViewExporter.hpp      ParaView输出
│                  NodalTxtExporter.hpp      文本输出
├── log/           Logger.hpp               日志封装
├── physics/       ElectrostaticFieldSolver.hpp
│                  ThermalFieldSolver.hpp
│                  MechanicalFieldSolver.hpp
├── post/          (后处理)
└── solver/        SolverFactory.hpp         求解器工厂

src/
├── app/           Application.cpp
├── assembly/      MfemPoissonAssembler.cpp
│                  MfemLinearElasticityAssembler.cpp
├── coeff/         CoefficientManager.cpp
├── coupling/      MultiPhysicsCoupler.cpp
├── frontend/      ConfigLoader.cpp
│                  Expression.cpp
├── integrator/    (自定义积分器实现)
├── io/            ParaViewExporter.cpp
│                  NodalTxtExporter.cpp
├── log/           Logger.cpp
├── physics/       ElectrostaticFieldSolver.cpp
│                  ThermalFieldSolver.cpp
│                  MechanicalFieldSolver.cpp
├── post/          (后处理实现)
└── solver/        SolverFactory.cpp
```

---

## 8. 验证结果

与 COMSOL 6.2 参考解对比：

### 8.1 Busbar 电→热 单向耦合（一阶单元）

- 配置文件：`busbar_electro_thermal_no_iteration.json`
- 网格：`busbar.msh`（一阶四面体）

| 指标    | 值        |
| ------- | --------- |
| MAE     | 3.75e-9 K |
| RMSE    | 4.28e-9 K |
| MaxAE   | 1.17e-8 K |
| MeanRel | 1.23e-11  |

**结论**：数值精度达到机器精度级别。

### 8.2 Busbar 电→热 单向耦合（二阶单元）

- 配置文件：`busbar_electro_thermal_no_iteration_order2.json`
- 网格：`busbar_order2.msh`（二阶四面体）

| 指标    | 值        |
| ------- | --------- |
| MAE     | 1.60e-6 K |
| MeanRel | 5.27e-9   |

### 8.3 Busbar 电↔热 Picard 迭代（温度相关电导率）

- 配置文件：`busbar_electro_thermal_iteration.json`
- 电导率表达式：$\sigma(T) = \frac{5.998 \times 10^7}{1 + 0.00393(T - 293.15)}$
- 收敛情况：4 次 Picard 迭代收敛（tolerance = 1e-6）

| 指标    | 值      |
| ------- | ------- |
| MAE     | 0.002 K |
| MeanRel | 6.68e-6 |

### 8.4 Bipolar Plate 热→力 耦合

- 配置文件：`bipolar_thermo_mechanical.json`
- 温度范围：[368.3, 439.9] K
- 位移分量对比：

| 分量  | MAE (m) | MeanRel |
| ----- | ------- | ------- |
| $u_x$ | 1.70e-6 | 2.6%    |
| $u_y$ | 1.19e-6 | 2.4%    |
| $u_z$ | 2.99e-7 | 17.2%*  |

\* $u_z$ 本身数值极小（~2e-6 m），绝对误差最低，相对误差受小分母影响偏大。

---

## 9. 扩展性考虑

### 9.1 添加瞬态求解

当前设计已预留扩展点：

1. `SimulationConfig` 可添加时间步进参数（`dt`, `t_end`, `time_scheme`）
2. 物理场求解器的 `Solve()` 接口无状态，可在外层时间循环中反复调用
3. `MultiPhysicsCoupler` 中的耦合循环可扩展为包裹时间步进
4. `Expression` 已支持 `t` 变量，时变边界条件无需改动

### 9.2 添加新物理场

1. 继承已有模式创建新的 `XxxFieldSolver`
2. 在 `Config.hpp` 添加对应的 `XxxFieldConfig`
3. 在 `ConfigLoader` 添加 JSON 解析逻辑
4. 在 `MultiPhysicsCoupler` 添加耦合策略

无需修改组装层或求解层代码。

### 9.3 性能进一步优化

- 对于大规模问题，可切换到 AMG 预条件的 CG 求解器
- 可引入 OpenMP 并行化高斯积分点上的系数求值
- 对于反复求解同一网格结构的问题，可缓存 `SparseMatrix` 的稀疏模式

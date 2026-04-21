# FEM Multi-Physics Solver

基于 [MFEM](https://mfem.org/) 构建的并行多物理场有限元求解器，支持静态/瞬态电-热-力耦合问题。

---

## 目录

1. [适用问题范围](#1-适用问题范围)
2. [数学模型与弱形式](#2-数学模型与弱形式)
3. [算法说明](#3-算法说明)
4. [项目结构层次](#4-项目结构层次)
5. [JSON 配置参考](#5-json-配置参考)
6. [API 说明](#6-api-说明)
7. [构建方法](#7-构建方法)
8. [运行示例](#8-运行示例)
9. [工具脚本](#9-工具脚本)

---

## 1. 适用问题范围

### 1.1 空间维度与单元类型

| 项目         | 支持情况                                              |
| ------------ | ----------------------------------------------------- |
| 空间维度     | 2D / 3D（由网格决定）                                 |
| 标量场单元   | H1-conforming Lagrange（阶次 `order` ≥ 1，默认 1 阶） |
| 向量场单元   | H1^d-conforming Lagrange（`order` 阶，用于力学场）    |
| 网格格式     | MFEM `.msh`（由 COMSOL `.mphtxt` 转换）、GMSH         |
| MPI 并行     | 所有求解器均并行；网格通过 METIS 自动分区             |
| 均匀网格细化 | 支持（`uniform_refinement_levels`）                   |

### 1.2 支持的物理场

#### 静电场（Electrostatic）
- 方程类型：稳态 Poisson 方程（Laplace 形式）
- 支持温度相关的电导率 σ(T)
- 边界条件：Dirichlet（固定电位）；时变 Dirichlet（`value_expr` 表达式，依赖变量 `t`）
- **不支持**：Neumann（直接电流注入）、Robin 边界

#### 热场（Thermal）
- 方程类型：稳态/瞬态热扩散方程（可含 Joule 热源）
- 支持温度相关的热导率 k(T)（通过表达式）
- 热容 ρCp 支持分段常数
- 边界条件：Dirichlet（固定温度）、Robin（对流换热 $h(T - T_\infty)$ 形式）
- 支持逐域热源（`source_by_domain`）和均匀热源（`source_default`）

#### 力学场（Mechanical）
- 方程类型：准静态线弹性（位移形式）
- 支持热应变（温度场驱动，单向耦合）
- 材料：各向同性线弹性（Young 模量 + Poisson 比 → Lamé 参数）
- 边界条件：
  - 位移 Dirichlet（三分量）
  - 法向位移约束（罚函数法）
  - 压力（Neumann）
  - 体力

### 1.3 支持的耦合模式

| 场组合                  | 耦合方式 | 备注                                    |
| ----------------------- | -------- | --------------------------------------- |
| E 单场                  | —        | 独立求解静电场                          |
| T 单场                  | —        | 独立求解热场（稳态或瞬态）              |
| M 单场                  | —        | 独立力学求解（无热场则热应变为零）      |
| E+T 单向（一次迭代）    | E→T      | 先解静电、再解热（不迭代）；适合弱耦合  |
| E+T 双向（Picard 迭代） | E⇌T      | σ(T) 反馈给 E；适合强耦合               |
| T+M 单向                | T→M      | 热场驱动力学热应变；热场先求解          |
| E+T+M                   | E⇌T → M  | 第一阶段 E-T 耦合，第二阶段力学批量求解 |

**Picard 判断逻辑**：当电导率表达式中包含 `T` 变量时，`ProjectConfig::NeedsPicardIteration()` 返回 `true`，启用双向迭代。

### 1.4 时间域范围

| 项目     | 支持情况                                                               |
| -------- | ---------------------------------------------------------------------- |
| 稳态     | 所有场组合                                                             |
| 瞬态     | 热场（BDF1/BDF2）+ 电场随时间变化的 BC；力学场在快照时刻批量准静态求解 |
| 时变 BC  | 电场 Dirichlet 支持 muParser 表达式（变量 `t`，单位秒）                |
| 初始条件 | 热场均匀初始温度 `initial_temperature`；力学初始位移为零               |

---

## 2. 数学模型与弱形式

### 2.1 静电场

**强形式**：
$$-\nabla \cdot (\sigma(T)\,\nabla V) = 0 \quad \text{in } \Omega$$

**弱形式**：
$$\int_\Omega \sigma(T)\,\nabla V \cdot \nabla \phi\,\mathrm{d}\Omega = 0 \quad \forall \phi \in H^1_0$$

### 2.2 热场（稳态）

**强形式**：
$$-\nabla \cdot (k\,\nabla T) = Q + \sigma|\nabla V|^2 \quad \text{in } \Omega$$
$$k\,\nabla T \cdot \mathbf{n} = q - h\,(T - T_\infty) \quad \text{on } \Gamma_R$$

**弱形式**：
$$\int_\Omega k\,\nabla T \cdot \nabla \psi\,\mathrm{d}\Omega + \int_{\Gamma_R} h\,T\,\psi\,\mathrm{d}S = \int_\Omega (Q + \sigma|\nabla V|^2)\,\psi\,\mathrm{d}\Omega + \int_{\Gamma_R} q\,\psi\,\mathrm{d}S$$

**矩阵形式**：$K T = F$，其中 $K = K_\mathrm{diff} + K_\mathrm{Robin}$

### 2.3 热场（瞬态）

**强形式**：
$$\rho C_p \frac{\partial T}{\partial t} - \nabla \cdot (k\,\nabla T) = Q$$

**矩阵半离散形式**：
$$C\,\dot{T} + K\,T = F(t)$$

- $K_{ij} = \int_\Omega k\,\nabla\phi_i \cdot \nabla\phi_j\,\mathrm{d}\Omega + \int_{\Gamma_R} h\,\phi_i\phi_j\,\mathrm{d}S$
- $C_{ij} = \int_\Omega \rho C_p\,\phi_i\phi_j\,\mathrm{d}\Omega$

### 2.4 力学场

**强形式**（准静态线弹性）：
$$-\nabla \cdot \boldsymbol{\sigma} = \mathbf{f} \quad \text{in } \Omega$$
$$\boldsymbol{\sigma} = \lambda\,(\mathrm{tr}\,\boldsymbol{\varepsilon})\mathbf{I} + 2\mu\,\boldsymbol{\varepsilon} - (3\lambda + 2\mu)\,\alpha_s(T - T_\mathrm{ref})\mathbf{I}$$

**弱形式**：
$$\int_\Omega \boldsymbol{\sigma} : \boldsymbol{\varepsilon}(\mathbf{v})\,\mathrm{d}\Omega = \int_\Omega \mathbf{f} \cdot \mathbf{v}\,\mathrm{d}\Omega + \int_{\Gamma_N} p\,\mathbf{v} \cdot \mathbf{n}\,\mathrm{d}S + \int_\Omega (3\lambda+2\mu)\,\alpha_s(T-T_\mathrm{ref})\,\nabla \cdot \mathbf{v}\,\mathrm{d}\Omega$$

法向位移约束通过罚函数处理：$K_\mathrm{stiff} \mathbf{u} = \mathbf{f}_\mathrm{body} + \mathbf{f}_\mathrm{pressure} + \mathbf{f}_\mathrm{thermal}$

---

## 3. 算法说明

### 3.1 线性代数求解器

通过 `fem::solver::CreateLinearSolver(name, ...)` 工厂函数统一创建，支持：

| 名称        | 类型                               | 实现类              | 适用场景                                   |
| ----------- | ---------------------------------- | ------------------- | ------------------------------------------ |
| `"pcg"`     | 迭代（并行 CG + BoomerAMG 预条件） | `MfemPcgSolver`     | 对称正定大系统（热、稳态电、力学）         |
| `"mumps"`   | 直接（并行 LU，MUMPS 库）          | `MfemMumpsSolver`   | 小中型系统、需要精确因子化（力学批量重用） |
| `"pardiso"` | 直接（Intel MKL PARDISO）          | `MfemPardisoSolver` | 仅限非对称矩阵（mtype=11）                 |
| `"umfpack"` | 直接（串行 UmfPack）               | `MfemUmfpackSolver` | 调试或串行场景                             |

**推荐组合**：
- 大型热场/静电场 → `"pcg"`
- 力学场（批量快照，LHS 不变）→ `"mumps"`（首次因子化后复用 `Mult()`）

### 3.2 瞬态时间积分

#### BDF1（Backward Euler）
$$\left(K + \frac{1}{\Delta t}C\right)T^{n+1} = F^{n+1} + \frac{1}{\Delta t}C\,T^n$$

系统矩阵 $A = K + \frac{1}{\Delta t}C$ 在 $\Delta t$ 不变时缓存复用（`EnsureSystemMatrix` 检查 `alpha_c` 是否改变）。

#### BDF2（变步长 BDF2）

步长比 $r = h_{n+1}/h_n$：
$$\alpha = \frac{1+2r}{(1+r)h},\quad \beta_n = \frac{1+r}{h},\quad \beta_{n-1} = \frac{r^2}{(1+r)h}$$

$$(K + \alpha C)T^{n+1} = F^{n+1} + C(\beta_n T^n - \beta_{n-1} T^{n-2})$$

退化为等步长时与标准 BDF2 $\left(\frac{3}{2h}\right)$ 一致。

### 3.3 自适应步长控制（预测-校正 WRMS）

#### 初始步长估计（Hairer-Wanner 风格）

冷启动时，计算 $\|C^{-1}(F - KT_0)\|_\infty$（对角线近似），得到 $\dot{T}_0$ 量级：
$$\Delta t_0 = \mathrm{clamp}\!\left(\frac{0.01\,d_0}{d_1},\; \Delta t_\mathrm{min},\; \Delta t_\mathrm{max}\right)$$
其中 $d_0 = \varepsilon_\mathrm{abs} + \varepsilon_\mathrm{rel}|T_0|$，$d_1 = \max(\|\dot{T}_0\|, d_0/\Delta t_\mathrm{max})$。

#### 误差估计

预测子（线性外推，history_depth ≥ 1）：
$$T_\mathrm{pred}^{n+1} = (1+r)\,T^n - r\,T^{n-1}, \quad r = \Delta t / \Delta t_\mathrm{prev}$$

校正子：BDF1（或 BDF2，若 history_depth ≥ 2）。

WRMS 误差：
$$E_\mathrm{WRMS} = \sqrt{\frac{1}{N_\mathrm{global}}\sum_{i=1}^N \left(\frac{T^\mathrm{corr}_i - T^\mathrm{pred}_i}{\varepsilon_\mathrm{abs} + \varepsilon_\mathrm{rel}|T^\mathrm{corr}_i|}\right)^2}$$

步长调整（$p$ = 积分阶次）：
$$\Delta t_\mathrm{new} = \Delta t \cdot \mathrm{clamp}\!\left(\eta_\mathrm{safe}\cdot E_\mathrm{WRMS}^{-1/(p+1)},\;\eta_\mathrm{min},\;\eta_\mathrm{max}\right)$$

$E_\mathrm{WRMS} \leq 1$ 时接受步，否则拒绝并将 $\Delta t$ 减半（或按公式缩小）。

#### BDF 阶次选择（order selection）

| `history_depth` | 预测子                           | 校正子      | 阶次选择                               |
| --------------- | -------------------------------- | ----------- | -------------------------------------- |
| 0（冷启动）     | 无                               | BDF1        | 总接受，步长扩大至 $\eta_\mathrm{max}$ |
| 1               | 线性（order-1）                  | BDF1        | 仅 order-1 WRMS                        |
| ≥ 2             | 二次（order-2）+ 线性（order-1） | BDF1 + BDF2 | 选 WRMS 更小的阶次和对应步长           |

二阶输出插值：使用三个连续快照的 Lagrange 二次多项式在输出时刻内插。

### 3.4 Picard 迭代（电热双向耦合）

```
for iter = 1, ..., max_iter:
    Solve E (使用当前 T 计算 σ(T))
    Solve T (使用当前 V 计算 Joule 热 σ|∇V|²)
    Δrel = ‖T_new - T_prev‖ / ‖T_prev‖
    if Δrel < tolerance: 收敛
    T_prev = relaxation * T_new + (1-relaxation) * T_prev
```

- **欠松弛**（`picard_relaxation < 1`）：适合 σ(T) 变化剧烈时
- **收敛准则**：true-DOF 空间下的 $L^2$ 相对变化量

### 3.5 力学场批量求解（瞬态第二阶段）

由于力学场为准静态（LHS $K$ 不随时间变化，仅 RHS 的热应变项随温度快照变化）：

1. **一次装配** `CacheStiffnessMatrix()`：组装刚度矩阵 $K$，仅执行一次
2. **首次求解** `SolveRHSOnly()`：调用 `CreateLinearSolver`，执行一次因子化
3. **后续快照** `SolveRHSOnly()`：调用 `solver->Mult(b, x)`（仅回代），避免重复因子化

直接求解器（MUMPS）在此模式下性价比最高。

### 3.6 矩阵缓存策略（热场）

- $K$（扩散矩阵）和 $C$（质量矩阵）在 `CacheMatrices()` 中各组装一次
- 系统矩阵 $A = \alpha_K K + \alpha_C C$ 仅在 $\alpha_C$ 改变时重建（自适应步长改变 $1/\Delta t$ 时触发）
- 稳态、BDF1、BDF2 各维护独立 `SolveState`，互不干扰

### 3.7 输出快照与插值

- 快照输出时刻由 `transient_output_interval` 均匀指定
- 计算步长与输出时刻一般不对齐：使用线性（2 点）或 Lagrange 二次（3 点）插值到输出时刻
- 输出内容：电位 $V$、温度 $T$、位移 $\mathbf{u}$（ParaView DataCollection + 文本文件）

---

## 4. 项目结构层次

```
FEM_mfem/
├── main.cpp                        # 入口：MPI 初始化 → Application::Run()
├── CMakeLists.txt
├── CMakePresets.json
│
├── include/fem/                    # 所有公共头文件
│   ├── app/
│   │   └── Application.hpp         # 顶层应用类
│   ├── frontend/
│   │   ├── Config.hpp              # 所有配置数据结构（ProjectConfig 等）
│   │   ├── ConfigLoader.hpp        # JSON → ProjectConfig 解析
│   │   └── Expression.hpp         # muParser 表达式包装（变量 x,y,z,t,T）
│   ├── physics/                    # 物理场求解器接口
│   │   ├── ElectrostaticFieldSolver.hpp
│   │   ├── ThermalFieldSolver.hpp
│   │   └── MechanicalFieldSolver.hpp
│   ├── coupling/
│   │   └── PicardCoupler.hpp       # E-T Picard 迭代器
│   ├── transient/
│   │   ├── TransientSolver.hpp     # 瞬态主控（BDF + 自适应 + 力学批量）
│   │   └── TransientUtils.hpp     # 快照辅助、WRMS 误差函数
│   ├── steady/
│   │   └── SteadySolver.hpp        # 稳态主控
│   ├── coeff/                      # 系数（材料属性 → MFEM Coefficient）
│   │   ├── CoefficientManager.hpp  # ExpressionCoefficient, JouleHeatingCoefficient 等
│   │   ├── ThermalCoeffs.hpp
│   │   ├── ElectrostaticCoeffs.hpp
│   │   └── MechanicalCoeffs.hpp
│   ├── assembler/                  # 矩阵/向量装配层
│   │   ├── ThermalAssembler.hpp
│   │   ├── ElectrostaticAssembler.hpp
│   │   └── MechanicalAssembler.hpp
│   ├── integrator/                 # 自定义 MFEM 积分器
│   │   ├── PoissonIntegrators.hpp  # 标量场：扩散、边界质量、体力/边界力
│   │   └── MechanicalIntegrators.hpp  # 热应变 LF、法向位移罚
│   ├── BCbuilder/                  # 边界条件构造
│   │   ├── ScalarBCBuilder.hpp     # 标量场 BC（Dirichlet + Robin）
│   │   └── MechanicalBCBuilder.hpp # 力学 BC（位移 + 法向约束 + 压力）
│   ├── solver/                     # 线性求解器抽象层
│   │   ├── ILinearSolver.hpp       # 接口（Solve + Mult）
│   │   ├── LinearSolverFactory.hpp # CreateLinearSolver 工厂
│   │   ├── MfemPcgSolver.hpp       # PCG + BoomerAMG
│   │   ├── MfemMumpsSolver.hpp     # MUMPS 直接求解器
│   │   ├── MfemPardisoSolver.hpp   # PARDISO 直接求解器
│   │   └── MfemUmfpackSolver.hpp   # UmfPack 直接求解器
│   ├── output/
│   │   ├── ResultExporter.hpp      # 并行结果导出（ParaView + txt）
│   │   ├── ParaviewExporter.hpp    # ParaView DataCollection
│   │   └── SolutionTextExporter.hpp
│   └── log/
│       └── Logger.hpp              # spdlog 封装（带时间戳）
│
├── src/                            # 实现文件，与 include/fem/ 一一对应
│   ├── app/Application.cpp
│   ├── frontend/{ConfigLoader,Expression}.cpp
│   ├── physics/{Electrostatic,Thermal,Mechanical}FieldSolver.cpp
│   ├── coupling/PicardCoupler.cpp
│   ├── transient/{TransientSolver,TransientUtils}.cpp
│   ├── steady/SteadySolver.cpp
│   ├── coeff/{CoefficientManager,Thermal,Electrostatic,Mechanical}Coeffs.cpp
│   ├── assembler/{Thermal,Electrostatic,Mechanical}Assembler.cpp
│   ├── integrator/{Poisson,Mechanical}Integrators.cpp
│   ├── BCbuilder/{Scalar,Mechanical}BCBuilder.cpp
│   ├── solver/{LinearSolverFactory,MfemPcg,MfemMumps,MfemPardiso,MfemUmfpack}Solver.cpp
│   ├── output/{ResultExporter,ParaviewExporter,SolutionTextExporter}.cpp
│   └── log/Logger.cpp
│
├── configs/                        # 示例配置文件
├── mesh_data/                      # 网格文件（.msh / .mphtxt）
├── results/                        # 输出目录（ParaView + txt）
├── docs/
│   └── physics.schema.json         # JSON Schema（IDE 补全用）
├── tools/
│   ├── mphtxt_to_msh.py           # COMSOL 网格转换
│   └── compare_solution_txt.py    # 与 COMSOL 参考解对比
└── mfem/                           # MFEM 源码（bundled 模式）
```

### 调用层次（稳态 E+T+M Picard 为例）

```
main()
└── Application::Run(config_path)
    ├── ConfigLoader::LoadFromFile()        → ProjectConfig
    └── SteadySolver::SolveSteady()
        ├── ElectrostaticFieldSolver(config)
        │   └── ElectrostaticCoeffs::BuildCoeffs()
        │       └── ExpressionCoefficient("electrical_conductivity", ...)
        ├── ThermalFieldSolver(config)
        │   └── ThermalCoeffs::BuildCoeffs()
        │       ├── ExpressionCoefficient("diffusion", ...)
        │       ├── PiecewiseConstantCoefficient("density"×"heat_capacity", ...)
        │       └── JouleHeatingCoefficient(&sigma_, &voltage_)
        ├── e_solver->SetTemperatureField(&T)   // 双指针晚绑定
        ├── t_solver->SetVoltageField(&V)
        ├── t_solver->CacheMatrices()           // 一次性组装 K, C
        ├── PicardCoupler::RunPicard(...)
        │   ├── e_solver->Solve()               // σ(T) 在 Eval 时实时取 T 值
        │   └── t_solver->SolveSteady()         // Joule 热 = σ|∇V|²
        ├── MechanicalFieldSolver(config)
        │   └── MechanicalCoeffs::BuildCoeffs()
        ├── m_solver->SetTemperatureField(&T)
        ├── m_solver->Solve()
        └── ResultExporter::ExportScalar/Vector(...)
```

### 晚绑定双指针模式（Late-bound Observer）

```cpp
// 内部持有指向"指针槽"的 const 指针
mfem::GridFunction * const *voltage_gf_;   // 读取 *voltage_gf_ 获得当前 GridFunction*

// 设置方式（外部把指针槽地址传入）
t_solver->SetVoltageField(&e_solver->GetVoltage());
// 此后 e_solver->GetVoltage() 更新时，Coefficient::Eval 自动读到最新值
```

该模式使 `JouleHeatingCoefficient` 等系数对象无需在每步重新构造，直接读取最新的 GridFunction 指针。

---

## 5. JSON 配置参考

所有配置通过单个 JSON 文件传入（`$schema` 字段可指向 `docs/physics.schema.json` 获得 IDE 补全）。

### 5.1 `simulation` 块

```jsonc
"simulation": {
  "mesh_path": "mesh_data/busbar.msh",      // 必填，相对于工作目录
  "order": 1,                               // 单元阶次，默认 1
  "uniform_refinement_levels": 0,           // 均匀细化次数，默认 0
  "output_dir": "results/my_case",          // 输出目录
  "log_level": "info",                      // trace/debug/info/warn/error

  // 线性求解器（可被各场 override）
  "solver": "pcg",                          // 全局默认：pcg | mumps | pardiso | umfpack
  "solver_electrostatic": "pcg",            // 仅覆盖静电场（可省略）
  "solver_thermal": "pcg",                  // 仅覆盖热场（可省略）
  "solver_mechanical": "mumps",             // 仅覆盖力学场（可省略）

  // Picard 迭代参数（仅 E-T 双向耦合时生效）
  "picard_max_iterations": 30,
  "picard_tolerance": 1e-6,                 // 相对 L2 变化量
  "picard_relaxation": 1.0,                 // 欠松弛系数 ∈ (0,1]

  // 瞬态控制
  "transient_enabled": true,
  "transient_t_start": 0.0,
  "transient_t_end": 0.01,
  "transient_dt": 0.001,                    // 固定步长（adaptive_dt=false 时使用）
  "transient_output_interval": 0.001,       // 输出快照间隔

  // 自适应步长控制（需 adaptive_dt=true）
  "adaptive_dt": true,
  "adaptive_reltol": 2.5e-6,               // WRMS 相对容差
  "adaptive_abstol": 5.0e-9,               // WRMS 绝对容差（K 量级）
  "dt_min": 1e-5,                          // 最小步长
  "dt_max": 0.003,                         // 最大步长（0=用 output_interval）
  "eta_safety": 0.9,                       // 安全因子
  "eta_max": 2.0,                          // 最大步长扩大倍数
  "eta_min": 0.2,                          // 最小步长缩小倍数

  // 可选：与 COMSOL 参考解对比
  "comsol_reference_path": "results/ref.txt",
  "compare_args": "--transient --fields-per-step 3 --field-names V,T,disp"
}
```

### 5.2 `materials` + `domain_materials` 块

```jsonc
"materials": [
  {
    "name": "copper",
    "properties": {
      // 标量场：热导率（可为纯数字或 muParser 表达式）
      "diffusion": "400.0",
      // 电导率：可依赖 T（变量名必须是 T，单位 K）
      "electrical_conductivity": "5.998e7/(1.0 + 0.0039*(T - 293.15))",
      // 力学场材料参数（单位：Pa）
      "young_modulus": "1.1e11",
      "poisson_ratio": "0.35",
      // 热膨胀系数（割线型，单位 1/K）
      "thermal_expansion_secant": "1.7e-5",
      // 热容（ρCp 直接提供，或分开提供 density + heat_capacity）
      "density": "8960",          // kg/m³
      "heat_capacity": "385"      // J/(kg·K)
    }
  }
],
"domain_materials": {
  "copper": [1],             // 域属性编号 1 对应 copper
  "titanium": [2, 3, 4]      // 域属性 2,3,4 对应 titanium
}
```

**支持的材料属性关键字**：

| 关键字                     | 用途                        | 是否支持表达式                    |
| -------------------------- | --------------------------- | --------------------------------- |
| `diffusion`                | 热导率 $k$                  | 是（可依赖 T,x,y,z）              |
| `electrical_conductivity`  | 电导率 $\sigma$             | 是（常依赖 T）                    |
| `density`                  | 密度 $\rho$                 | 否（常数）                        |
| `heat_capacity`            | 比热容 $C_p$                | 否（常数，与 density 相乘得 ρCp） |
| `young_modulus`            | 弹性模量 $E$                | 否                                |
| `poisson_ratio`            | 泊松比 $\nu$                | 否                                |
| `thermal_expansion_secant` | 割线热膨胀系数 $\alpha_s$   | 否                                |
| `lambda_lame`              | Lamé 第一参数（可替代 E+ν） | 否                                |
| `mu_lame`                  | 剪切模量（可替代 E+ν）      | 否                                |

> 注：若同时提供 `young_modulus`+`poisson_ratio`，程序在 `BuildCoeffs` 中自动换算为 Lamé 参数。若直接提供 `lambda_lame`+`mu_lame`，则跳过换算。

### 5.3 `fields` 数组

#### 静电场字段

```jsonc
{
  "type": "electrostatic",
  "source_default": "0.0",           // 体电荷源（通常为 0）
  "reference_temperature": 293.15,   // σ(T) 的参考温度（单场模式或初始 Picard 温度）

  "dirichlet_boundary_conditions": [
    {
      "bdr_attributes": [43],
      "value": 0.02                  // 固定值（V）
    },
    {
      "bdr_attributes": [43],
      "value_expr": "0.02*sin(100*pi*t)"  // 时变（优先级高于 value）
    },
    {
      "bdr_attributes": [8, 15],
      "value": 0.0                   // 接地
    }
  ]
}
```

#### 热场字段

```jsonc
{
  "type": "thermal",
  "source_default": "0.0",             // 均匀热源（W/m³）
  "source_by_domain": {
    "3": "1.1086e7"                    // 域 3 的热源（覆盖 source_default）
  },
  "initial_temperature": 293.15,       // 初始/参考温度（K）

  "dirichlet_boundary_conditions": [
    { "bdr_attributes": [1], "value": 300.0 }
  ],

  "robin_boundary_conditions": [
    {
      "bdr_attributes": [2, 3, 4],
      // BC 形式：k ∂T/∂n = q - l*T  （即 h*(T∞ - T)，其中 l=h, q=h*T∞）
      "l": 5.0,                        // 对流系数 h（W/(m²·K)）
      "q": 1465.75                     // h * T∞（W/m²）
    }
  ]
}
```

> Robin BC 说明：程序将 $h(T_\infty - T)$ 分解为 $q - l \cdot T$，故 `l` = $h$，`q` = $h \cdot T_\infty$。

#### 力学场字段

```jsonc
{
  "type": "mechanical",
  "reference_temperature": 293.15,    // 热应变参考温度（K）

  "body_force_default": [0.0, 0.0, 0.0],  // 均匀体力（N/m³）

  "displacement_boundary_conditions": [
    {
      "bdr_attributes": [1, 2],
      "value": [0.0, 0.0, 0.0]        // 全约束（可分量分别施加）
    }
  ],

  "normal_displacement_boundary_conditions": [
    {
      "bdr_attributes": [3, 5, 20],
      "penalty": 1.0e16               // 罚参数（越大越近似于法向零位移）
    }
  ],

  "pressure_boundary_conditions": [
    {
      "bdr_attributes": [24],
      "value": 1.0e6                  // 压力（Pa，正值为压缩方向）
    }
  ]
}
```

---

## 6. API 说明

### 6.1 顶层入口 `fem::app::Application`

```cpp
// include/fem/app/Application.hpp
class Application {
public:
    int Run(const std::string &config_path);
};
```

**说明**：唯一公共入口。读取 JSON、调度稳态或瞬态求解器、执行 COMSOL 对比。

**典型使用**（`main.cpp`）：
```cpp
MPI_Init(&argc, &argv);
fem::log::Init("info");
fem::app::Application app;
int ret = app.Run(argv[1]);
MPI_Finalize();
```

---

### 6.2 配置加载 `fem::frontend::ConfigLoader`

```cpp
// include/fem/frontend/ConfigLoader.hpp
static ProjectConfig LoadFromFile(const std::string &path);
```

**说明**：解析 JSON → 填充 `ProjectConfig`（含材料数据库、FEConfig 中的 mesh/FESpace）。

---

### 6.3 物理场求解器

#### `fem::physics::ElectrostaticFieldSolver`

```cpp
explicit ElectrostaticFieldSolver(frontend::ProjectConfig &config);

// 耦合接口：在 Picard 循环前调用一次
void SetTemperatureField(mfem::GridFunction *temperature_gf);

// 可选：注入外部已构建的 sigma（通常由 ThermalFieldSolver 共享）
void SetElectricalConductivity(mfem::Coefficient *sigma);

// 更新时变 Dirichlet BC（值通过 value_expr 重新求值）
void UpdateBoundaryConditions(double t);

// 求解（自动组装、施加 BC、求解线性系统）
void Solve();

// 获取结果
mfem::ParGridFunction &GetVoltage();

// 获取当前 sigma 系数（可传给 ThermalFieldSolver）
mfem::Coefficient *GetSigma();
```

**推荐使用模式**（E+T Picard）：
```cpp
auto e = std::make_unique<ElectrostaticFieldSolver>(config);
auto t = std::make_unique<ThermalFieldSolver>(config);
e->SetTemperatureField(&t->GetTemperature());  // 晚绑定 T
t->SetVoltageField(&e->GetVoltage());          // 晚绑定 V
t->SetElectricalConductivity(e->GetSigma());   // 共享 sigma
```

---

#### `fem::physics::ThermalFieldSolver`

```cpp
explicit ThermalFieldSolver(frontend::ProjectConfig &config);
// 初始化：温度设为 initial_temperature，构建 BC 数据

// 耦合接口（耦合时调用，单场可跳过）
void SetVoltageField(mfem::GridFunction *voltage_gf);
void SetElectricalConductivity(mfem::Coefficient *sigma);

// 矩阵缓存（必须在任意 Solve* 前调用一次）
void CacheMatrices();

// 稳态求解
void SolveSteady();

// 瞬态：设置 dt 和 T^n（每步调用）
void EnableTransient(double dt, mfem::GridFunction *T_old);

// BDF1（Backward Euler）
void SolveBDF1();

// BDF2（变步长，需提供 T^{n-1} 和 dt_{prev}）
void SolveBDF2(const mfem::Vector &T_nm1_true, double dt_prev);

// 手动组装当前热源 RHS（F_current_），通常由 Solve* 内部调用
void AssembleSourceRHS();

// 计算 max |dT/dt|_0 用于自适应初始步长估计
double ComputeInitialRate();  // 需先 CacheMatrices()

// 获取结果
mfem::ParGridFunction &GetTemperature();
const mfem::Vector &GetBDF2Solution() const;  // BDF2 解的 true-DOF 向量
```

**瞬态使用模式**：
```cpp
t_solver->CacheMatrices();
// 时间步循环：
T_old = t_solver->GetTemperature();
t_solver->EnableTransient(dt, &T_old);
t_solver->SolveBDF1();         // 或 SolveBDF2(T_nm1, dt_prev)
```

---

#### `fem::physics::MechanicalFieldSolver`

```cpp
explicit MechanicalFieldSolver(frontend::ProjectConfig &config);

// 热耦合（可选，不设置则热应变为零）
void SetTemperatureField(mfem::GridFunction *temperature_gf);

// 完整求解（内部调用 CacheStiffness + SolveRHSOnly）
void Solve();

// 批量模式：先缓存刚度矩阵（一次）
void CacheStiffnessMatrix();

// 批量模式：仅组装 RHS 并求解（每个快照调用，复用因子化）
void SolveRHSOnly();

// 获取结果
mfem::ParGridFunction &GetDisplacement();
```

**批量快照推荐模式**：
```cpp
m_solver->SetTemperatureField(&t_solver->GetTemperature());
m_solver->CacheStiffnessMatrix();         // 因子化一次
for (auto &snap : snapshots) {
    // 将快照温度设置到 temperature GridFunction
    t_solver->GetTemperature().Distribute(snap.temperature);
    m_solver->SolveRHSOnly();             // 仅回代
}
```

---

### 6.4 耦合控制 `fem::coupling::PicardCoupler`

```cpp
// include/fem/coupling/PicardCoupler.hpp
explicit PicardCoupler(const frontend::SimulationConfig &sim_config);

PicardResult RunPicard(
    physics::ElectrostaticFieldSolver *e_solver,   // 可为 nullptr
    physics::ThermalFieldSolver *t_solver,         // 可为 nullptr
    bool has_e, bool has_t,
    PicardSolveMode mode   // Steady 或 BackwardEuler
);

struct PicardResult { int iterations; bool converged; };
```

**模式**：
- `Steady`：调用 `t_solver->SolveSteady()`
- `BackwardEuler`：调用 `t_solver->SolveBDF1()`（需先 `EnableTransient()`）

**单场降级**：若 `has_e=false`，直接调用 `t_solver->SolveSteady/BDF1()`，返回 `{1, true}`；反之亦然。

---

### 6.5 瞬态主控 `fem::transient::TransientSolver`

```cpp
explicit TransientSolver(frontend::ProjectConfig &config);
void SolveTransient();  // 运行完整瞬态模拟
```

内部自动执行：
1. 创建三个场求解器
2. 设置耦合接口（晚绑定指针）
3. 一致初始步长估计（`ComputeInitialRate()`）
4. 自适应 BDF1/BDF2 时间推进 + Picard 每步耦合
5. 线性/二次插值生成输出快照
6. 力学场批量准静态求解（阶段二）
7. 调用 `ResultExporter::ExportTransient()`

---

### 6.6 系数管理 `fem::coeff`

#### `ExpressionCoefficient`

```cpp
// 基于 muParser 的分段表达式系数
ExpressionCoefficient(
    const std::string &property_name,  // 如 "electrical_conductivity"
    const frontend::MaterialDatabase &db,
    mfem::Mesh &mesh,
    mfem::GridFunction * const *temperature_gf = nullptr  // 可选晚绑定温度
);
void SetReferenceTemperature(double T_ref);  // 无温度场时的备用值
```

**表达式变量**：`x`, `y`, `z`（积分点坐标，m），`T`（温度，K），`t`（时间，s）。

#### `JouleHeatingCoefficient`

```cpp
JouleHeatingCoefficient(
    mfem::Coefficient * const *sigma = nullptr,
    mfem::GridFunction * const *voltage = nullptr
);
// Eval 时计算 σ * |∇V|² at integration point
```

#### `CombinedHeatSourceCoefficient`

将体热源 $Q_\mathrm{base}$ 和 Joule 热 $Q_J$ 叠加：

```cpp
CombinedHeatSourceCoefficient(mfem::Coefficient *base_source,
                               JouleHeatingCoefficient *joule_source);
```

---

### 6.7 线性求解器 `fem::solver`

```cpp
// 工厂函数
std::unique_ptr<ILinearSolver> CreateLinearSolver(
    const std::string &solver_name,   // "pcg" | "mumps" | "pardiso" | "umfpack"
    double rel_tol, double abs_tol,
    int max_iter, int print_level
);

// 接口
class ILinearSolver {
    virtual void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) = 0;
    // 直接求解器：复用上次因子化；迭代求解器：重新求解
    virtual void Mult(const mfem::Vector &b, mfem::Vector &x) = 0;
};
```

---

### 6.8 自适应误差工具 `fem::transient`

```cpp
// 计算 WRMS 误差并给出建议步长
WRMSResult ComputeAdaptiveError(
    const mfem::Vector &T_pred,   // 预测子
    const mfem::Vector &T_corr,   // 校正子
    double dt,
    double reltol, double abstol,
    double eta_safety, double eta_max, double eta_min,
    double dt_min,
    int order = 1   // BDF 阶次
);

struct WRMSResult {
    double wrms;     // WRMS 误差（≤1 为接受）
    double dt_new;   // 建议下一步长
    bool accept;     // 是否接受本步
};
```

---

### 6.9 日志 `fem::log`

```cpp
void Init(const std::string &level);          // "trace"/"debug"/"info"/"warn"/"error"
std::shared_ptr<spdlog::logger> Get();         // 获取全局 logger

// 用法（任意 .cpp 中）
auto logger = fem::log::Get();
logger->info("Phase 1: step={}, dt={:.4e}", step, dt);
logger->debug("T=[{:.4e}, {:.4e}]", T.Min(), T.Max());
```

输出格式带有相对时间戳（秒），便于性能分析。

---

## 7. 构建方法

### 依赖

| 依赖                         | 版本要求 | 获取方式                                      |
| ---------------------------- | -------- | --------------------------------------------- |
| MFEM（含 MPI、METIS、HYPRE） | ≥ 4.6    | 内置 `mfem/` 子目录（bundled 模式）或系统安装 |
| OpenMPI                      | ≥ 4.0    | 系统包管理器                                  |
| MUMPS                        | 推荐     | MFEM 编译时开启                               |
| spdlog                       | ≥ 1.9    | FetchContent 自动拉取                         |
| nlohmann_json                | ≥ 3.10   | FetchContent 自动拉取                         |
| muparser                     | ≥ 2.3    | FetchContent 自动拉取                         |
| CMake                        | ≥ 3.16   | —                                             |
| Ninja                        | 推荐     | —                                             |

### 构建步骤

```bash
# Debug 构建（内置 MFEM）
cmake --preset mfem-bundled-linux-debug
cmake --build build/linux-debug -j$(nproc) --target FEM_solution

# Release 构建
cmake --preset linux-release
cmake --build build/linux-release -j$(nproc) --target FEM_solution
```

可用预设名称（`CMakePresets.json` 中定义）：
- `mfem-bundled-linux-debug` — Debug，使用 `mfem/` 源码
- `mfem-system-linux-debug` — Debug，使用系统已安装 MFEM
- `linux-release` — Release `-O3`

---

## 8. 运行示例

### 8.1 稳态电热 Picard（强耦合）

```bash
mpirun -np 8 ./build/linux-release/bin/FEM_solution \
    ./configs/busbar_electro_thermal_picard.json
```

输出：`results/busbar_picard/`

### 8.2 瞬态电热力（自适应步长 + 时变 BC）

```bash
mpirun -np 8 ./build/linux-release/bin/FEM_solution \
    ./configs/busbar_etm_nonlinear_transient_vbc.json
```

配置特点：
- 电导率 $\sigma(T)$ 温度依赖（Picard 耦合）
- 电位 BC 为正弦时变：`"0.2*sin(100*pi*t)"`
- 自适应步长：`adaptive_reltol=2.5e-6`，`adaptive_abstol=5e-9`
- 力学场批量准静态求解（使用 MUMPS）

### 8.3 热-力耦合（稳态，含压力边界）

```bash
mpirun -np 8 ./build/linux-release/bin/FEM_solution \
    ./configs/bipolar_thermo_mechanical.json
```

### 8.4 典型输出

```
[0.000s] [info] ==================================
[0.000s] [info]   FEM Multi-Physics Solver
[0.000s] [info] ==================================
[0.012s] [info] === Transient solve: t=[0, 0.01], dt=0.001, interval=0.001 ===
[0.012s] [info]   Adaptive PC: reltol=2.50e-06 abstol=5.00e-09 dt=[1.00e-05,3.00e-03]
[0.015s] [info] Phase 1: Electro-thermal time stepping
[0.020s] [info] Consistent init: max |dT/dt|_0=1.3241e+02 K/s, dt_init=2.4500e-04
[1.736s] [info] Phase 1 done: 19 accepted, 5 rejected, 11 snapshots
[1.736s] [info] Phase 2: Mechanical batch solve (11 snapshots)
[3.945s] [info] Phase 2 done: 11 mechanical solves
[3.945s] [info] Transient complete: 19 steps, 5 rejected, 11 snapshots
[5.489s] [info] All done.
```

---

## 9. 工具脚本

### 9.1 网格转换 `tools/mphtxt_to_msh.py`

```bash
python3 tools/mphtxt_to_msh.py mesh_data/busbar.mphtxt mesh_data/busbar.msh
```

将 COMSOL `.mphtxt` 导出格式转换为 MFEM/GMSH 兼容的 `.msh` 格式。

### 9.2 结果对比 `tools/compare_solution_txt.py`

```bash
python3 tools/compare_solution_txt.py \
    --mine results/busbar_transient_vbc/transient.txt \
    --comsol results/busbar_transient_vbc.txt \
    --transient --fields-per-step 3 --field-names V,T,disp
```

输出各时间步各场的最大绝对误差与相对误差。可通过配置中的 `comsol_reference_path` 和 `compare_args` 字段自动触发。

---

## 附录：材料属性名称对照

| JSON 属性名                | 物理量         | 符号       | 单位     |
| -------------------------- | -------------- | ---------- | -------- |
| `diffusion`                | 热导率         | $k$        | W/(m·K)  |
| `electrical_conductivity`  | 电导率         | $\sigma$   | S/m      |
| `density`                  | 密度           | $\rho$     | kg/m³    |
| `heat_capacity`            | 比热容         | $C_p$      | J/(kg·K) |
| `young_modulus`            | 弹性模量       | $E$        | Pa       |
| `poisson_ratio`            | 泊松比         | $\nu$      | —        |
| `lambda_lame`              | Lamé 第一参数  | $\lambda$  | Pa       |
| `mu_lame`                  | 剪切模量       | $\mu$      | Pa       |
| `thermal_expansion_secant` | 割线热膨胀系数 | $\alpha_s$ | 1/K      |

**Lamé 参数换算**：
$$\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \qquad \mu = \frac{E}{2(1+\nu)}$$

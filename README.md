# FEM Multi-Physics Solver

基于 MFEM 的 MPI 并行多物理场有限元求解器。当前项目主线是显式矩阵装配 + Hypre 并行线性代数，支持静电、热传导和准静态线弹性，以及电-热、热-力、电-热-力组合计算。

本文档按当前代码实现整理。项目接口、MFEM 常用接口、MFEM MPI 接口和自适应时间步说明分别见 [docs/project_interfaces.md](docs/project_interfaces.md)、[docs/mfem_common_interfaces.md](docs/mfem_common_interfaces.md)、[docs/mfem_mpi_interfaces.md](docs/mfem_mpi_interfaces.md) 与 [docs/adaptive_timestepping.md](docs/adaptive_timestepping.md)；原有完整索引 [docs/mfem_api.md](docs/mfem_api.md) 暂时保留。

## 目录

1. [当前能力](#1-当前能力)
2. [物理模型与弱形式](#2-物理模型与弱形式)
3. [求解流程](#3-求解流程)
4. [JSON 配置](#4-json-配置)
5. [构建](#5-构建)
6. [运行示例](#6-运行示例)
7. [项目结构](#7-项目结构)
8. [工具脚本](#8-工具脚本)
9. [当前边界](#9-当前边界)

## 1. 当前能力

### 1.1 离散与并行

| 项目 | 当前实现 |
| ---- | -------- |
| 空间维度 | 2D / 3D，由网格决定 |
| 标量场空间 | `mfem::H1_FECollection(order, dim)` + `mfem::ParFiniteElementSpace` |
| 矢量场空间 | `H1^dim`，用于力学位移 |
| 网格 | MFEM/Gmsh `.msh`；仓库中有 COMSOL `.mphtxt` 转换脚本 |
| MPI | 始终创建 `ParMesh` 和并行有限元空间，`mpirun -np 1` 也是并行代码路径 |
| 线性系统 | `HypreParMatrix` + true-DOF 向量 |
| 输出 | ParaView DataCollection 与 rank 0 聚合后的节点文本文件 |

`ConfigLoader::BuildFE()` 先读取串行网格并做可选均匀细化，然后用 `mfem::ParMesh(MPI_COMM_WORLD, mesh)` 建立分布式网格。项目已经没有单独的串行求解主路径。

### 1.2 物理场组合

| 组合 | 求解方式 |
| ---- | -------- |
| E | 独立静电场 |
| T | 独立热场 |
| M | 独立准静态力学场，未设置温度场时热应变为零 |
| E + T | 若电导率不依赖 `T`，先 E 后 T；若依赖 `T`，启用 Picard 迭代 |
| T + M | 先 T 后 M，温度驱动热应变 |
| E + T + M | 先电热耦合，再力学求解；瞬态中力学在输出快照上批量求解 |

Picard 是否启用由 `ProjectConfig::NeedsPicardIteration()` 决定：只要任一材料的 `electrical_conductivity` 表达式实际使用变量 `T`，电热场就按双向耦合处理。

### 1.3 稳态与瞬态

| 项目 | 当前实现 |
| ---- | -------- |
| 稳态 | E、T、M 及其组合 |
| 瞬态 | 热场时间推进为主，电场可通过时变 Dirichlet BC 随时间更新 |
| 时间积分 | 后向欧拉 BDF1；自适应步长时可在历史足够后尝试变步长 BDF2 |
| 误差控制 | 预测-校正 WRMS 误差范数 |
| 输出时刻 | 按 `transient_output_interval` 生成，当前使用相邻快照线性插值 |
| 力学瞬态 | 准静态批量求解：缓存刚度矩阵，每个输出快照只重建 RHS |

## 2. 物理模型与弱形式

### 2.1 静电场

强形式：

$$
-\nabla \cdot \left(\sigma(T)\nabla V\right) = 0
\quad \text{in } \Omega
$$

Robin 边界约定：

$$
\sigma \nabla V \cdot \mathbf{n} = q - lV
\quad \text{on } \Gamma_R
$$

弱形式：

$$
\int_{\Omega} \sigma(T)\nabla V \cdot \nabla \phi\,d\Omega+ \int_{\Gamma_R} lV\phi\,dS
= \int_{\Gamma_R} q\phi\,dS
$$

当前静电场支持：

- 分域电导率 `electrical_conductivity`，可为 muparser 表达式。
- 表达式变量：`x`, `y`, `z`, `T`。
- Dirichlet 边界，`value` 可为数字或依赖 `t` 的字符串表达式。
- Robin 边界，形式与标量通用装配一致：`sigma dV/dn = q - l V`。

注意：静电场的 `source_default` 和 `source_by_domain` 当前没有接入系数构造，体源实际固定为 0。

### 2.2 热场

稳态强形式：

$$
-\nabla \cdot (k\nabla T) =
Q_\mathrm{base} + \sigma |\nabla V|^2
\quad \text{in } \Omega
$$

Robin 边界：

$$
k\nabla T \cdot \mathbf{n} = q - lT
\quad \text{on } \Gamma_R
$$

稳态弱形式：

$$
\int_{\Omega} k\nabla T \cdot \nabla \psi\,d\Omega+ \int_{\Gamma_R} lT\psi\,dS= \int_{\Omega} (Q_\mathrm{base}+\sigma|\nabla V|^2)\psi\,d\Omega+ \int_{\Gamma_R} q\psi\,dS
$$

瞬态强形式：

$$
\rho C_p\frac{\partial T}{\partial t}- \nabla \cdot (k\nabla T)= Q_\mathrm{base} + \sigma|\nabla V|^2
$$

瞬态半离散形式：

$$
C\dot{T} + KT = F(t)
$$

其中：

$$
K_{ij} =
\int_{\Omega} k\nabla\phi_i\cdot\nabla\phi_j\,d\Omega+ \int_{\Gamma_R} l\phi_i\phi_j\,dS,
\qquad
C_{ij} = \int_{\Omega}\rho C_p\phi_i\phi_j\,d\Omega
$$

当前热场支持：

- 热导率 `diffusion`：当前按分域常数读取。
- 热容 `rho_cp`：可直接提供，也可由 `density * heat_capacity` 自动生成。
- 基础热源 `source_default` / `source_by_domain`：当前按常数字符串读取。
- Joule 热源：若启用静电场，自动叠加 `sigma |grad V|^2`。
- Dirichlet 与 Robin 边界。

若要表达对流换热 `h (T_inf - T)`，配置中使用 `l = h`，`q = h * T_inf`。

### 2.3 力学场

准静态线弹性强形式：

$$
-\nabla\cdot\boldsymbol{\sigma} = \mathbf{f}
\quad \text{in } \Omega
$$

$$
\boldsymbol{\sigma}= \lambda\,\mathrm{tr}(\boldsymbol{\varepsilon})\mathbf{I}+ 2\mu\,\boldsymbol{\varepsilon}- (3\lambda+2\mu)\alpha_s(T-T_\mathrm{ref})\mathbf{I}
$$

弱形式按“弹性刚度 + 体力 + 压力边界 + 热应变载荷”拆分：

$$
\int_\Omega
\left(\lambda\,\mathrm{tr}(\varepsilon(\mathbf{u}))\mathbf{I}+ 2\mu\,\varepsilon(\mathbf{u})\right)
:\varepsilon(\mathbf{v})\,d\Omega=
\int_\Omega \mathbf{f}\cdot\mathbf{v}\,d\Omega+ \int_{\Gamma_p} p\,\mathbf{v}\cdot\mathbf{n}\,dS+ \int_\Omega
(3\lambda+2\mu)\alpha_s(T-T_\mathrm{ref})
\nabla\cdot\mathbf{v}\,d\Omega
$$

法向位移弱约束通过罚函数加入刚度项：

$$
\gamma\int_{\Gamma_n}
(\mathbf{u}\cdot\mathbf{n})(\mathbf{v}\cdot\mathbf{n})\,dS
$$

当前力学场支持：

- 各向同性线弹性。
- `young_modulus` + `poisson_ratio` 自动换算为 `lambda_lame` 和 `mu_lame`。
- 位移 Dirichlet。
- 法向位移罚约束 `normal_displacement_boundary_conditions`。
- 压力边界 `pressure_boundary_conditions`。
- 分域或全局体力。
- 温度场驱动的热应变；未设置温度场时使用参考温度，热应变为零。

## 3. 求解流程

### 3.1 线性求解器

配置项 `solver` 是全局默认值，`solver_electrostatic`、`solver_thermal`、`solver_mechanical` 可逐场覆盖。

| 名称 | 实现 | 当前路径 |
| ---- | ---- | -------- |
| `pcg` | `mfem::CGSolver` + `mfem::HypreBoomerAMG` | 并行迭代求解，要求 `HypreParMatrix` |
| `mumps` | `mfem::MUMPSSolver` | 并行直接求解，当前设置为 SPD 矩阵类型 |
| `pardiso` | `mfem::CPardisoSolver` | MKL Cluster PARDISO，并行直接求解 |

仓库中保留了 `MfemUmfpackSolver` 包装类，但工厂函数当前不暴露 `umfpack`，且该实现会抛出“串行 UMFPACK 不支持并行配置”的异常。配置文件中请使用 `pcg`、`mumps` 或 `pardiso`。

### 3.2 稳态流程

稳态由 `SteadySolver::SolveSteady()` 编排：

```text
Load config and FE spaces
  -> solve electro-thermal stage
     -> E only, T only, E->T one-way, or E<->T Picard
  -> solve mechanical stage if enabled
     -> optionally read final temperature for thermal strain
  -> export scalar/vector results
```

电热 Picard 每轮执行：

```text
solve E with current T
solve T with current V and sigma
check true-DOF relative L2 change of T
apply relaxation if needed
```

### 3.3 瞬态流程

瞬态由 `TransientSolver::SolveTransient()` 编排：

```text
create solvers
set E/T late-bound coupling pointers
cache thermal K and C
solve initial E if present
estimate initial dt when adaptive_dt=true
while output times remain:
    set T^n and current dt
    update time-varying electrostatic Dirichlet values
    solve backward-Euler Picard step
    estimate WRMS error when adaptive_dt=true
    optionally solve BDF2 and choose accepted order
    accept/reject step and update history
    linearly interpolate output snapshots
batch mechanical solves over output snapshots
export transient results
```

BDF1 system:

```text
(K + C / dt) T^{n+1} = F^{n+1} + (C / dt) T^n
```

Variable-step BDF2 system, with `r = dt / dt_prev`:

```text
alpha = (1 + 2r) / ((1 + r) dt)
beta_n = (1 + r) / dt
beta_nm1 = r^2 / ((1 + r) dt)

(K + alpha C) T^{n+1}
  = F^{n+1} + C (beta_n T^n - beta_nm1 T^{n-1})
```

预测-校正误差使用 true-DOF 空间上的 WRMS 范数：

$$
E_\mathrm{WRMS}=
\sqrt{
\frac{1}{N_\mathrm{global}}
\sum_i
\left(
\frac{|T_i^\mathrm{corr}-T_i^\mathrm{pred}|}
{\varepsilon_\mathrm{abs}
+\varepsilon_\mathrm{rel}|T_i^\mathrm{corr}|}
\right)^2
}
$$

步长调整公式：

$$
\Delta t_\mathrm{new}=
\Delta t\cdot
\mathrm{clamp}
\left(
\eta_\mathrm{safety}E_\mathrm{WRMS}^{-1/(p+1)},
\eta_\mathrm{min},
\eta_\mathrm{max}
\right)
$$

阶次选择逻辑是效率优先：如果 BDF1 和 BDF2 都满足 WRMS 接受准则，选择建议下一步长更大的那个阶次。

### 3.4 缓存策略

- 热场在 `CacheMatrices()` 中组装并缓存 `K` 与 `C`。
- 热场每个求解模式维护独立 `SolveState`，系统矩阵仅在 `alpha_c` 改变时重建。
- 力学场在批量模式下缓存 stiffness `ParBilinearForm` 和求解器。后续快照仅重新装配 RHS，并调用 `solver->Mult()` 复用直接求解器的分解或迭代器状态。
- 边界 marker 与 essential true dofs 在 solver 构造时缓存。

## 4. JSON 配置

所有配置由单个 JSON 文件传入。示例位于 [configs](configs)，schema 位于 [docs/physics.schema.json](docs/physics.schema.json)。

### 4.1 `simulation`

```jsonc
{
  "simulation": {
    "mesh_path": "mesh_data/busbar.msh",
    "order": 1,
    "uniform_refinement_levels": 0,
    "output_dir": "results/my_case",
    "log_level": "info",

    "solver": "pcg",
    "solver_electrostatic": "pcg",
    "solver_thermal": "pcg",
    "solver_mechanical": "mumps",

    "picard_max_iterations": 30,
    "picard_tolerance": 1e-6,
    "picard_relaxation": 1.0,

    "transient_enabled": true,
    "transient_t_start": 0.0,
    "transient_t_end": 0.01,
    "transient_dt": 0.001,
    "transient_output_interval": 0.001,

    "adaptive_dt": true,
    "adaptive_reltol": 2.5e-6,
    "adaptive_abstol": 5.0e-9,
    "dt_min": 1e-5,
    "dt_max": 0.003,
    "eta_safety": 0.9,
    "eta_max": 2.0,
    "eta_min": 0.2,

    "comsol_reference_path": "results/reference.txt",
    "compare_args": "--transient --fields-per-step 3 --field-names V,T,disp",
    "curve_reference_path": "results/reference.txt",
    "curve_compare_args": "--random-points 8 --fields-per-step 3 --field-names V,T,disp --plot-dir results/curve_plots"
  }
}
```

`dt_max = 0` 时，代码使用 `transient_output_interval` 作为有效最大步长。瞬态启用时，加载器会检查自适应步长容差与步长上下限。

### 4.2 `materials` 与 `domain_materials`

```jsonc
{
  "materials": [
    {
      "name": "copper",
      "properties": {
        "diffusion": "400.0",
        "electrical_conductivity": "5.998e7/(1.0 + 0.0039*(T - 293.15))",
        "density": "8960",
        "heat_capacity": "385",
        "young_modulus": "1.1e11",
        "poisson_ratio": "0.35",
        "thermal_expansion_secant": "1.7e-5"
      }
    }
  ],
  "domain_materials": {
    "copper": [1],
    "titanium": [2, 3, 4, 5, 6, 7]
  }
}
```

| 属性 | 用途 | 当前读取方式 |
| ---- | ---- | ------------ |
| `diffusion` | 热导率 `k` | 分域常数 |
| `rho_cp` | 体积热容 | 分域常数 |
| `density`, `heat_capacity` | 自动生成 `rho_cp` | 常数，使用 `std::stod` |
| `electrical_conductivity` | 电导率 `sigma` | muparser 表达式，可依赖 `x,y,z,T` |
| `young_modulus`, `poisson_ratio` | 自动生成 Lamé 参数 | 常数 |
| `lambda_lame`, `mu_lame` | 直接提供 Lamé 参数 | 常数 |
| `thermal_expansion_secant` | 割线热膨胀系数 | 常数 |

### 4.3 静电场

```jsonc
{
  "type": "electrostatic",
  "dirichlet_boundary_conditions": [
    { "bdr_attributes": [43], "value": "0.2*sin(100*pi*t)" },
    { "bdr_attributes": [8, 15], "value": 0.0 }
  ],
  "robin_boundary_conditions": [
    { "bdr_attributes": [10], "l": 1.0, "q": 0.0 }
  ],
  "reference_temperature": 293.15
}
```

时变边界值使用 `value` 字段传字符串，不是单独的 `value_expr` 字段。字符串表达式会在瞬态步开始时按当前时间 `t` 更新；空间变量此处不会传入真实边界坐标。

### 4.4 热场

```jsonc
{
  "type": "thermal",
  "source_default": "0.0",
  "source_by_domain": {
    "3": "1.1086e7"
  },
  "initial_temperature": 293.15,
  "dirichlet_boundary_conditions": [
    { "bdr_attributes": [1], "value": 300.0 }
  ],
  "robin_boundary_conditions": [
    { "bdr_attributes": [2, 3, 4], "l": 5.0, "q": 1465.75 }
  ]
}
```

当热场存在时，`initial_temperature` 会传播为静电场与力学场的参考温度。

### 4.5 力学场

```jsonc
{
  "type": "mechanical",
  "reference_temperature": 293.15,
  "body_force_default": [0.0, 0.0, 0.0],
  "displacement_boundary_conditions": [
    { "bdr_attributes": [8, 15, 43], "value": [0.0, 0.0, 0.0] }
  ],
  "normal_displacement_boundary_conditions": [
    { "bdr_attributes": [3, 5, 20], "penalty": 1.0e16 }
  ],
  "pressure_boundary_conditions": [
    { "bdr_attributes": [24], "value": 1.0e6 }
  ]
}
```

## 5. 构建

### 5.1 主要依赖

| 依赖 | 说明 |
| ---- | ---- |
| CMake | `CMakeLists.txt` 要求 3.16；`CMakePresets.json` 要求 3.23 |
| MPI | Linux preset 使用 `mpicc` / `mpicxx` |
| MFEM | 默认使用仓库内 `mfem/`，也可切换系统安装版 |
| HYPRE / METIS | MFEM MPI 路径需要 |
| MUMPS / MKL CPARDISO | 直接求解器可选但当前 CMake 默认尝试开启 |
| spdlog / nlohmann_json / muparser | 优先 `find_package`，找不到时通过 `FetchContent` 获取 |

### 5.2 Debug 构建

```bash
cmake --preset mfem-bundled-linux-debug
cmake --build build/linux-debug --target FEM_solution -j$(nproc)
```

### 5.3 Release 构建

```bash
cmake --preset linux-release
cmake --build build/linux-release --target FEM_solution -j$(nproc)
```

也可以直接构建 CMake 目标名：

```bash
cmake --build build/linux-release --target my_fem -j$(nproc)
```

最终可执行文件名是 `FEM_solution`，位于对应 build 目录的 `bin/` 下。

## 6. 运行示例

### 6.1 稳态电热单向耦合

```bash
mpirun -np 4 ./build/linux-release/bin/FEM_solution \
  configs/busbar_electro_thermal_oneway.json
```

### 6.2 稳态电热 Picard

```bash
mpirun -np 4 ./build/linux-release/bin/FEM_solution \
  configs/busbar_electro_thermal_picard.json
```

### 6.3 瞬态电热力，自适应步长与时变电压

```bash
mpirun -np 4 ./build/linux-release/bin/FEM_solution \
  configs/busbar_etm_nonlinear_transient_vbc.json
```

该配置包含：

- `electrical_conductivity` 对 `T` 的依赖，触发电热 Picard。
- 电压边界 `"0.2*sin(100*pi*t)"`。
- 自适应 WRMS 步长控制。
- `solver_mechanical = "mumps"`，用于力学批量快照求解。

### 6.4 稳态热-力耦合

```bash
mpirun -np 4 ./build/linux-release/bin/FEM_solution \
  configs/bipolar_thermo_mechanical.json
```

### 6.5 输出文件

稳态输出通常包含：

```text
results/<case>/
  electrostatic_*.pvtu / thermal_*.pvtu / mechanical_*.pvtu
  voltage.txt
  temperature.txt
  displacement.txt
```

瞬态输出包含每个场的 ParaView collection，以及宽表格式：

```text
results/<case>/transient.txt
```

若配置了 `comsol_reference_path`，求解结束后 rank 0 会自动调用 [tools/compare_solution_txt.py](tools/compare_solution_txt.py) 做坐标匹配误差对比。若配置了 `curve_reference_path`，还会调用 [tools/compare_time_curves.py](tools/compare_time_curves.py) 对采样点瞬态时间曲线做对比。

## 7. 项目结构

```text
FEM_mfem/
├── main.cpp
├── CMakeLists.txt
├── CMakePresets.json
├── include/fem/
│   ├── app/              # Application 顶层入口
│   ├── frontend/         # Config, ConfigLoader, Expression
│   ├── physics/          # E/T/M 场求解器
│   ├── coupling/         # PicardCoupler
│   ├── transient/        # TransientSolver 与工具函数
│   ├── steady/           # SteadySolver
│   ├── coeff/            # MFEM Coefficient 封装
│   ├── assembler/        # 三个物理场的装配层
│   ├── integrator/       # 自定义积分子
│   ├── BCbuilder/        # 边界条件 marker 与系数构造
│   ├── solver/           # 线性求解器抽象与工厂
│   ├── output/           # ParaView 与 txt 导出
│   └── log/              # spdlog 封装
├── src/                  # 与 include/fem/ 对应的实现
├── configs/              # 示例配置
├── docs/                 # schema、接口文档、MFEM API、自适应步长、并行分析文档
├── tools/                # 网格转换与结果对比脚本
├── mesh_data/            # 本地网格数据，已被 .gitignore 忽略
└── mfem/                 # bundled MFEM 源码
```

主调用链：

```text
main
  -> mfem::Mpi::Init, mfem::Hypre::Init
  -> Application::Run
     -> ConfigLoader::LoadFromFile
        -> BuildFE: Mesh -> ParMesh -> ParFiniteElementSpace
     -> SteadySolver or TransientSolver
     -> ResultExporter
     -> optional COMSOL comparison
```

## 8. 工具脚本

### 8.1 COMSOL 网格转换

```bash
python3 tools/mphtxt_to_msh.py mesh_data/busbar.mphtxt mesh_data/busbar.msh
```

脚本读取 COMSOL `.mphtxt` 网格并输出 Gmsh/MFEM 可读的 `.msh`。二阶三角形、二阶四面体等节点顺序在脚本内有专门处理。

### 8.2 结果对比

```bash
python3 tools/compare_solution_txt.py \
  --mine results/busbar_transient_vbc/transient.txt \
  --comsol results/busbar_transient_vbc.txt \
  --transient --fields-per-step 3 --field-names V,T,disp
```

脚本按坐标容差匹配节点，输出 MAE、RMSE、MaxAE、MeanRel 等指标。

采样点时间曲线对比：

```bash
python3 tools/compare_time_curves.py \
  --mine results/busbar_transient_vbc/transient.txt \
  --comsol results/busbar_transient_vbc.txt \
  --random-points 8 --seed 2026 \
  --fields-per-step 3 --field-names V,T,disp \
  --plot-dir results/busbar_transient_vbc/curve_plots
```

也可以用 `--points "x,y,z; x,y,z"` 或 `--point-file points.txt` 指定采样点。`--plot-dir` 会为每个采样点和字段输出一张 FEM/COMSOL 曲线对比 SVG；`--output-csv` 可同时导出逐时间步误差明细。

## 9. 当前边界

- 热场 `diffusion`、`rho_cp` 与基础热源当前按常数读取；温度相关热导率尚未接入热场系数构造。
- 静电场材料电导率支持表达式，常用变量是 `T`；时变电压边界表达式通过 `value` 字符串提供，当前只传入时间 `t`。
- 当前瞬态主路径面向含热场的问题；已有示例均为 E+T+M。若要做纯热瞬态或 E+T 瞬态，需要先补齐导出阶段对缺失 vector FE space 的保护。
- `umfpack` 包装类保留但不属于当前可配置求解器。
- 当前没有 MPI 子通信器分组、场间任务并发、GPU/Partial Assembly 主路径，也没有 rank 内自定义 OpenMP 装配。

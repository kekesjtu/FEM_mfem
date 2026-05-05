# 项目接口文档（不含第三方库接口）

本文档只整理本项目自身暴露和约定的接口：命令行入口、JSON 配置、模块职责、输出文件和工具脚本输入输出。底层数值库、并行库、日志库、JSON 库和表达式库的类型与 API 不在本文展开；相关细节分别见 [mfem_api.md](mfem_api.md)、[mfem_common_interfaces.md](mfem_common_interfaces.md) 与 [mfem_mpi_interfaces.md](mfem_mpi_interfaces.md)。

## 1. 入口接口

### 1.1 命令行

可执行文件接受一个可选配置路径：

```bash
FEM_solution [config_path]
```

| 参数 | 必填 | 说明 |
| ---- | ---- | ---- |
| `config_path` | 否 | JSON 配置文件路径。若省略，程序使用源码中写定的默认路径。当前仓库建议显式传入 `configs/*.json`。 |

推荐运行方式：

```bash
mpirun -np 4 ./build/linux-release/bin/FEM_solution configs/busbar_etm_nonlinear_transient_vbc.json
```

### 1.2 应用流程

顶层流程由应用入口统一编排：

```text
读取 JSON 配置
  -> 初始化日志
  -> 构造网格、空间和材料/边界配置
  -> 根据 transient_enabled 选择稳态或瞬态流程
  -> 导出结果
  -> 可选调用结果对比脚本
```

应用层只关心“配置 -> 求解 -> 输出”的生命周期，不直接暴露底层矩阵、向量或并行对象。

## 2. JSON 配置接口

配置文件是项目最主要的用户接口。顶层结构如下：

```jsonc
{
  "simulation": {},
  "materials": [],
  "domain_materials": {},
  "fields": []
}
```

| 顶层字段 | 必填 | 说明 |
| -------- | ---- | ---- |
| `simulation` | 是 | 网格、阶次、输出、求解器、耦合迭代、瞬态和自动对比参数。 |
| `materials` | 否 | 材料属性表。 |
| `domain_materials` | 否 | 材料名到网格域编号列表的映射。 |
| `fields` | 是 | 启用的物理场列表。 |

完整 schema 见 [physics.schema.json](physics.schema.json)。

### 2.1 `simulation`

| 字段 | 类型 | 默认值 | 说明 |
| ---- | ---- | ------ | ---- |
| `mesh_path` | string | 无 | 网格文件路径。 |
| `order` | integer | `1` | 有限元阶次。 |
| `uniform_refinement_levels` | integer | `0` | 求解前均匀细化次数。 |
| `output_dir` | string | `results` | 输出目录。 |
| `log_level` | string | `info` | 日志等级。 |
| `solver` | string | `pcg` | 全局默认线性求解器。 |
| `solver_electrostatic` | string | `solver` | 静电场求解器覆盖项。 |
| `solver_thermal` | string | `solver` | 热场求解器覆盖项。 |
| `solver_mechanical` | string | `solver` | 力学场求解器覆盖项。 |
| `picard_max_iterations` | integer | `30` | 非线性电热 Picard 最大迭代数。 |
| `picard_tolerance` | number | `1e-6` | Picard 收敛阈值。 |
| `picard_relaxation` | number | `1.0` | Picard 温度松弛因子。 |
| `transient_enabled` | boolean | `false` | 是否启用瞬态流程。 |
| `transient_t_start` | number | `0.0` | 起始时间。 |
| `transient_t_end` | number | `1.0` | 结束时间。 |
| `transient_dt` | number | `1.0` | 初始或固定时间步长。 |
| `transient_output_interval` | number | `10.0` | 瞬态输出快照间隔。 |
| `adaptive_dt` | boolean | `false` | 是否启用自适应时间步。 |
| `adaptive_reltol` | number | `1e-3` | 自适应误差相对容差。 |
| `adaptive_abstol` | number | `1e-6` | 自适应误差绝对容差。 |
| `dt_min` | number | `1e-6` | 自适应最小步长。 |
| `dt_max` | number | `0.0` | 自适应最大步长；`0` 表示使用输出间隔。 |
| `eta_safety` | number | `0.9` | 步长缩放安全系数。 |
| `eta_max` | number | `2.0` | 单次最大放大倍数。 |
| `eta_min` | number | `0.2` | 单次最小缩小倍数。 |
| `comsol_reference_path` | string | `""` | 求解后自动全场对比的参考文件；为空则跳过。 |
| `compare_args` | string | `""` | 传给全场对比脚本的额外参数。 |
| `curve_reference_path` | string | `""` | 求解后自动采样点时间曲线对比的参考文件；为空则跳过。 |
| `curve_compare_args` | string | `""` | 传给时间曲线对比脚本的额外参数。 |

### 2.2 材料与域映射

材料属性以字符串存储，便于表达常数或表达式：

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
    "copper": [1, 2]
  }
}
```

| 属性 | 用途 |
| ---- | ---- |
| `diffusion` | 热导率。 |
| `rho_cp` | 体积热容。若缺省，可由 `density * heat_capacity` 自动生成。 |
| `electrical_conductivity` | 电导率，可依赖 `x,y,z,T`。 |
| `young_modulus`, `poisson_ratio` | 自动生成 Lamé 参数。 |
| `lambda_lame`, `mu_lame` | 直接提供 Lamé 参数。 |
| `thermal_expansion_secant` | 热膨胀系数。 |

### 2.3 物理场配置

`fields` 是数组，按 `type` 区分物理场。

| `type` | 说明 |
| ------ | ---- |
| `electrostatic` | 静电场。支持 Dirichlet、Robin、电导率表达式和时变电压边界。 |
| `thermal` | 热场。支持基础热源、Joule 热源、Dirichlet、Robin 和初始温度。 |
| `mechanical` | 准静态线弹性。支持位移边界、法向位移罚约束、压力边界和体力。 |

标量场通用边界配置：

```jsonc
{
  "dirichlet_boundary_conditions": [
    { "bdr_attributes": [1], "value": 293.15 }
  ],
  "robin_boundary_conditions": [
    { "bdr_attributes": [2, 3], "l": 5.0, "q": 1465.75 }
  ]
}
```

力学场边界配置：

```jsonc
{
  "displacement_boundary_conditions": [
    { "bdr_attributes": [1], "value": [0.0, 0.0, 0.0] }
  ],
  "normal_displacement_boundary_conditions": [
    { "bdr_attributes": [2], "penalty": 1.0e16 }
  ],
  "pressure_boundary_conditions": [
    { "bdr_attributes": [3], "value": 1.0e6 }
  ]
}
```

## 3. 模块接口

### 3.1 前端配置层

| 模块 | 输入 | 输出 | 职责 |
| ---- | ---- | ---- | ---- |
| `Application` | 配置路径 | 进程返回码 | 顶层生命周期编排。 |
| `ConfigLoader` | JSON 文件 | 项目配置对象 | 解析配置、材料、物理场与离散设置，并做基本校验。 |
| `Expression` | 表达式字符串 | 数值结果 | 计算依赖 `x,y,z,t,T` 的用户表达式，并检测变量依赖。 |

配置层不负责物理装配，只负责把用户输入变成内部可消费的数据结构。

### 3.2 物理场层

| 模块 | 输入 | 输出 | 职责 |
| ---- | ---- | ---- | ---- |
| `ElectrostaticFieldSolver` | 项目配置、可选温度场 | 电势场 | 组装并求解静电方程，处理时变电压边界。 |
| `ThermalFieldSolver` | 项目配置、可选电势场和电导率 | 温度场 | 求解稳态热方程或 BDF 时间步，支持矩阵缓存。 |
| `MechanicalFieldSolver` | 项目配置、可选温度场 | 位移场 | 求解准静态线弹性问题，支持批量快照的刚度矩阵缓存。 |
| `PicardCoupler` | 电场求解器、热场求解器、耦合参数 | 迭代次数与收敛状态 | 在电导率依赖温度时执行电热 Picard 迭代。 |

### 3.3 求解流程层

| 模块 | 输入 | 输出 | 职责 |
| ---- | ---- | ---- | ---- |
| `SteadySolver` | 项目配置 | 稳态结果文件 | 编排 E/T/M 及其组合的稳态求解。 |
| `TransientSolver` | 项目配置 | 瞬态结果文件 | 编排热主导的时间推进、快照缓存、插值输出和瞬态导出。 |
| `TransientUtils` | 相邻快照、误差参数 | 输出快照或步长建议 | 保存场快照、插值输出时刻、计算自适应误差。 |

### 3.4 边界、系数和输出层

| 模块 | 输入 | 输出 | 职责 |
| ---- | ---- | ---- | ---- |
| `ScalarBCBuilder` | 标量场边界配置 | 标量边界数据 | 构造 Dirichlet/Robin 边界标记与系数数据。 |
| `MechanicalBCBuilder` | 力学边界配置 | 力学边界数据 | 构造位移、法向位移罚约束和压力边界数据。 |
| `CoefficientManager` | 材料库、域映射、场变量 | 系数对象 | 统一材料属性、表达式系数和 Joule 热源。 |
| `ResultExporter` | 输出目录、场结果、时间信息 | ParaView 与文本文件 | 导出稳态和瞬态结果。 |
| `SolutionTextExporter` | 场结果、网格节点、时间序列 | 文本宽表 | 生成与 COMSOL 对齐的节点文本格式。 |

## 4. 输出文件接口

### 4.1 稳态文本

稳态标量场：

```text
% ...
% x y z value
x y z value
```

稳态矢量场：

```text
% ...
% x y z ux uy uz magnitude
x y z ux uy uz magnitude
```

常见文件名：

| 文件 | 内容 |
| ---- | ---- |
| `voltage.txt` | 电势节点值。 |
| `temperature.txt` | 温度节点值。 |
| `displacement.txt` | 位移分量与模长。 |

### 4.2 瞬态文本宽表

瞬态结果写入：

```text
results/<case>/transient.txt
```

格式：

```text
% x y z V @ t=t0 T @ t=t0 disp @ t=t0 V @ t=t1 T @ t=t1 disp @ t=t1 ...
x y z V0 T0 disp0 V1 T1 disp1 ...
```

约定：

- 前三列固定为 `x y z`。
- 后续列按时间步排列。
- 每个时间步默认有三个场：`V,T,disp`。
- 输出快照时间来自 `transient_output_interval`，自适应步长时可由相邻已接受快照插值得到。

## 5. 工具脚本接口

### 5.1 网格转换

```bash
python3 tools/mphtxt_to_msh.py input.mphtxt output.msh
```

输入 COMSOL 文本网格，输出 Gmsh/MFEM 可读网格。

### 5.2 全场结果对比

```bash
python3 tools/compare_solution_txt.py --mine result.txt --comsol reference.txt
```

常用参数：

| 参数 | 说明 |
| ---- | ---- |
| `--mine` | 本项目导出的文本结果。 |
| `--comsol` | COMSOL 导出的参考文本结果。 |
| `--value-col` | 强制比较指定值列。 |
| `--transient` | 按瞬态宽表解析。 |
| `--fields-per-step` | 每个时间步的字段数。 |
| `--field-names` | 字段名列表。 |
| `--tol` | 坐标匹配容差。 |

### 5.3 采样点时间曲线对比

```bash
python3 tools/compare_time_curves.py --mine transient.txt --comsol reference.txt --random-points 8
```

常用参数：

| 参数 | 说明 |
| ---- | ---- |
| `--random-points N` | 从可匹配节点中随机选取 `N` 个点。 |
| `--points "x,y,z;..."` | 显式指定采样点坐标。 |
| `--point-file file.txt` | 从文件读取采样点坐标。 |
| `--fields-per-step` | 每个时间步字段数，默认 `3`。 |
| `--field-names` | 字段名列表，默认 `V,T,disp`。 |
| `--fields` | 只比较部分字段，可填字段名或索引。 |
| `--output-csv` | 导出逐点、逐时刻误差明细。 |
| `--plot-dir` | 输出 FEM 与 COMSOL 时间曲线对比图片。 |
| `--plot-format` | 图片格式，默认 `svg`；`png` 需要 matplotlib。 |
| `--tol` | 坐标匹配容差。 |

自动调用时，配置 `curve_reference_path` 和 `curve_compare_args` 即可。

## 6. 错误与返回约定

| 场景 | 行为 |
| ---- | ---- |
| 配置文件打不开或 JSON 解析失败 | 抛出异常，主程序记录 fatal 日志并返回非零码。 |
| 缺少必需配置段 | 抛出异常。 |
| 自适应时间步参数非法 | 抛出异常。 |
| 自动对比参考文件不存在 | 记录 warning，跳过该对比，不中断求解结果保存。 |
| Python 对比脚本返回非零码 | 记录 warning，主求解流程已经完成。 |

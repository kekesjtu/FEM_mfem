# MFEM API 参考文档（本项目完整版）

> 本文档详尽记录本项目 (`FEM_mfem`) 使用的所有 MFEM 接口及其数据流，
> 特别聚焦于 **VSize / TrueVSize / GlobalTrueVSize** 三层 DOF 体系，
> 以及 MPI 多进程架构在有限元求解中的协作机制。

---

## 目录

1. [MPI 多进程架构总览](#1-mpi-多进程架构总览)
2. [DOF 三层体系：VSize → TrueVSize → GlobalTrueVSize](#2-dof-三层体系)
3. [网格与几何](#3-网格与几何)
4. [有限元空间](#4-有限元空间)
5. [GridFunction — 场变量载体](#5-gridfunction--场变量载体)
6. [双线性形式与线性形式](#6-双线性形式与线性形式)
7. [矩阵与向量](#7-矩阵与向量)
8. [积分子（Integrators）](#8-积分子integrators)
9. [系数（Coefficients）](#9-系数coefficients)
10. [线性求解器](#10-线性求解器)
11. [边界条件处理](#11-边界条件处理)
12. [完整求解数据流](#12-完整求解数据流)
13. [输出与导出](#13-输出与导出)
14. [MPI 通讯点汇总](#14-mpi-通讯点汇总)

---

## 1. MPI 多进程架构总览

### 1.1 什么是 MPI 多进程

MPI（Message Passing Interface）并行不是多线程——每个 MPI rank 是一个**独立的操作系统进程**，拥有自己的内存空间。它们之间看不到彼此的变量，必须通过**显式消息传递**来交换数据。

一个生动的比喻：

> 想象一个巨大的拼图（有限元网格），你找了 4 个朋友来拼。
> 你把拼图切成 4 块（METIS 分区），每人拿走一块到各自的房间（独立进程）。
> 每个人只能看到自己那块拼图（本地网格分区）。
> 当需要让边界处的图案对齐时，他们必须通过对讲机（MPI 通讯）互相报告边界颜色。
> 最终每人拼好自己的部分，合在一起就是完整答案。

### 1.2 本项目的 MPI 启动流程

```
main.cpp
├── mfem::Mpi::Init(argc, argv)    // 初始化 MPI 运行时
├── mfem::Hypre::Init()            // 初始化 Hypre 并行线性代数库
└── Application::Run()
    └── ConfigLoader::BuildFE()
        ├── Mesh("busbar.msh")     // rank 0 加载完整网格
        └── ParMesh(MPI_COMM_WORLD, *mesh)  // METIS 自动分区，每个 rank 持有一份子网格
```

运行命令示例：
```bash
# 单进程
mpirun -np 1 ./build/linux-release/bin/FEM_solution configs/busbar_electro_thermal_oneway.json

# 并行（2 个进程）
mpirun -np 2 ./build/linux-release/bin/FEM_solution configs/busbar_electro_thermal_oneway.json
```

### 1.3 MPI 架构示意图

```
┌────────────────────────────────────────────────────────────────┐
│                     MPI_COMM_WORLD                             │
│                                                                │
│  ┌──────────────────┐        ┌──────────────────┐              │
│  │    Rank 0         │        │    Rank 1         │              │
│  │                  │        │                  │              │
│  │  ParMesh (分区0) │        │  ParMesh (分区1) │              │
│  │  ┌────────────┐  │        │  ┌────────────┐  │              │
│  │  │ 本地单元    │  │        │  │ 本地单元    │  │              │
│  │  │ + 鬼魂层    │  │        │  │ + 鬼魂层    │  │              │
│  │  └────────────┘  │        │  └────────────┘  │              │
│  │                  │        │                  │              │
│  │  ParFESpace      │        │  ParFESpace      │              │
│  │  VSize = 120     │        │  VSize = 115     │              │
│  │  TrueVSize = 100 │        │  TrueVSize = 95  │              │
│  │                  │        │                  │              │
│  │  GridFunction    │        │  GridFunction    │              │
│  │  (size = VSize)  │        │  (size = VSize)  │              │
│  │                  │        │                  │              │
│  │  HypreParMatrix  │        │  HypreParMatrix  │              │
│  │  (行 = TrueVSize)│        │  (行 = TrueVSize)│              │
│  └──────────────────┘        └──────────────────┘              │
│           ↕ MPI_Allreduce / ParallelAssemble ↕                 │
│                                                                │
│              GlobalTrueVSize = 100 + 95 = 195                  │
│              (全局真实自由度 = 各 rank TrueVSize 之和)           │
└────────────────────────────────────────────────────────────────┘
```

### 1.4 每个 rank 持有什么

| 数据             | 每个 rank 独立持有的内容                      | 跨 rank 共享/通讯                  |
| ---------------- | --------------------------------------------- | ---------------------------------- |
| `ParMesh`        | 本地子网格 + 鬼魂层 (ghost layer)             | METIS 分区后自动构建               |
| `ParFESpace`     | 本地 DOF 编号 → 全局 DOF 映射                 | 构建时自动协调共享节点             |
| `GridFunction`   | 本地 VSize 长度的向量（含共享节点的冗余副本） | `GetTrueDofs/SetFromTrueDofs` 转换 |
| `HypreParMatrix` | 本地行对应 TrueVSize 行的 CSR 块              | Hypre 内部管理列的全局索引         |
| `Vector` (RHS)   | TrueVSize 长度                                | `ParLinearForm::ParallelAssemble`  |

---

## 2. DOF 三层体系

这是理解 MFEM 并行最关键的概念。三个 "Size" 代表三种不同的自由度视角：

### 2.1 三层定义

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│  VSize (L-vector 大小)                                          │
│  ═══════════════════                                            │
│  "本地所有 DOF，包含共享节点的冗余副本"                           │
│                                                                 │
│  • 串行: VSize = 网格上所有节点 × 每节点自由度数                  │
│  • 并行: VSize = 本地节点 + 鬼魂节点(共享边界冗余) × 每节点DOF数  │
│  • GridFunction 的大小 = VSize                                  │
│  • 用于：物理空间插值、后处理、可视化                            │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  TrueVSize (T-vector 大小)                                      │
│  ═══════════════════════                                        │
│  "本进程独自拥有的 DOF，去掉共享节点的冗余"                       │
│                                                                 │
│  • 串行: TrueVSize = VSize（没有冗余，二者相等）                  │
│  • 并行: TrueVSize ≤ VSize（共享节点只算在其中一个 rank 上）      │
│  • HypreParMatrix 的行数 = TrueVSize                            │
│  • 用于：线性系统求解、矩阵-向量运算                            │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  GlobalTrueVSize (全局大小)                                      │
│  ═══════════════════════════                                     │
│  "整个问题的真实自由度总数 = 所有 rank 的 TrueVSize 之和"         │
│                                                                 │
│  • 各 rank 的 TrueVSize 互不重叠，直接求和即可                    │
│  • 仅用于日志输出，不在计算中直接使用                             │
│  • 等价于串行情况下的全局 VSize                                   │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 生动比喻

> 想象一栋公寓楼的信箱系统：
>
> - **VSize** = 每层楼看到的信箱数。如果你住在两层楼的交界处（共享节点），
>   两层楼都给你留了一个信箱，所以"楼层看到的信箱数"比实际住户多。
>
> - **TrueVSize** = 你这一层楼真正负责投递的信箱数。
>   共享节点的信箱只由其中一层楼负责，另一层虽然有冗余信箱但不负责投递。
>
> - **GlobalTrueVSize** = 整栋楼实际住户总数。
>   把每层楼负责的信箱数加起来就是了，因为职责不重叠。

### 2.3 Prolongation 和 Restriction — 两个空间之间的桥梁

```
                    Prolongation (P)
  T-vector ──────────────────────────────→ L-vector
  (TrueVSize)    "展开到含冗余副本"         (VSize)
                                              │
                     Restriction (R = Pᵀ)     │
  T-vector ←──────────────────────────────── L-vector
  (TrueVSize)    "收集并求和冗余副本"         (VSize)
```

**本项目实际代码**（[ThermalFieldSolver.cpp](src/physics/ThermalFieldSolver.cpp) `ApplyBCAndSolve`）：

```cpp
// ① 从 L-vector (GridFunction) 提取 T-vector
mfem::Vector x_true(N);    // N = TrueVSize
temperature_.GetTrueDofs(x_true);  // L → T (restriction)

// ② 在 T-vector 空间中求解线性系统
state.solver->Solve(A, b, x_true);  // A 是 HypreParMatrix, 行数 = TrueVSize

// ③ 将解从 T-vector 展开回 L-vector
output->SetSize(par_fespace->GetVSize());  // VSize
par_fespace->GetProlongationMatrix()->Mult(x_true, *output);  // T → L (prolongation)
```

### 2.4 本项目中的具体数值示例

以 `busbar.msh`（6 面体棱柱网格）+ 2 MPI ranks 为例：

```
串行模式:
  Mesh: elements=1260, vertices=1660
  Scalar FE space (order 1): VSize = TrueVSize = 1660
  Vector FE space (order 1): VSize = TrueVSize = 4980 (1660 × 3)

并行模式 (np=2):
  Rank 0: local elements ≈ 630
    Scalar: VSize ≈ 910,  TrueVSize ≈ 830
  Rank 1: local elements ≈ 630
    Scalar: VSize ≈ 920,  TrueVSize ≈ 830
  
  GlobalTrueVSize = 830 + 830 = 1660  ← 和串行一致！
  
  差异说明:
    VSize(910) > TrueVSize(830) 因为分区边界上约有 80 个共享节点，
    这些节点在两个 rank 上各有一份冗余副本。
```

---

## 3. 网格与几何

### 3.1 Mesh — 串行网格

| 接口                           | 说明                       | 本项目使用位置                                                                                                     |
| ------------------------------ | -------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `Mesh(filename, gen, ref)`     | 从文件加载网格             | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp) `BuildFE()`                                                      |
| `UniformRefinement()`          | 均匀细化                   | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp)                                                                  |
| `Dimension()`                  | 获取空间维度 (2D/3D)       | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp)                                                                  |
| `GetNE()`                      | 单元总数                   | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp) 日志输出                                                         |
| `GetNV()`                      | 顶点总数                   | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp), [SolutionTextExporter.cpp](src/output/SolutionTextExporter.cpp) |
| `GetNBE()`                     | 边界单元数                 | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp)                                                                  |
| `bdr_attributes.Max()`         | 最大边界属性编号           | 所有 FieldSolver 中构建 BC marker                                                                                  |
| `GetElementVertices(e, verts)` | 获取第 e 个单元的顶点列表  | [SolutionTextExporter.cpp](src/output/SolutionTextExporter.cpp)                                                    |
| `GetElementTransformation(e)`  | 获取参考→物理坐标变换      | [SolutionTextExporter.cpp](src/output/SolutionTextExporter.cpp)                                                    |
| `GetElementBaseGeometry(e)`    | 单元几何类型 (tet/hex/etc) | [SolutionTextExporter.cpp](src/output/SolutionTextExporter.cpp)                                                    |

### 3.2 ParMesh — 并行网格

| 接口                         | 说明                                                 | 本项目使用位置                                         |
| ---------------------------- | ---------------------------------------------------- | ------------------------------------------------------ |
| `ParMesh(comm, serial_mesh)` | 从串行网格 + MPI 通讯域创建并行网格，自动 METIS 分区 | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp)      |
| `GetComm()`                  | 获取 MPI 通讯域                                      | [ResultExporter.cpp](src/output/ResultExporter.cpp)    |
| `GetNE()`                    | 本分区的单元数                                       | [ConfigLoader.cpp](src/frontend/ConfigLoader.cpp) 日志 |

`ParMesh` 继承自 `Mesh`，所有 `Mesh` 接口均可使用。构造时 METIS 自动将全局网格分为 `WorldSize` 块，每个 rank 只持有一块。

### 3.3 几何工具

| 接口                                 | 说明                         |
| ------------------------------------ | ---------------------------- |
| `mfem::Geometries.GetVertices(geom)` | 获取参考单元顶点对应的积分点 |
| `IntegrationPoint`                   | 参考空间中的一个点 (ξ, η, ζ) |
| `IntegrationRule`                    | 一组积分点 + 权重            |
| `ElementTransformation`              | 参考空间 ↔ 物理空间的变换    |

---

## 4. 有限元空间

### 4.1 H1_FECollection

```cpp
auto fec = std::make_unique<mfem::H1_FECollection>(order, dim);
```

H1 连续有限元族（拉格朗日型），所有场（标量/矢量）均使用此族。

| 参数    | 含义                         |
| ------- | ---------------------------- |
| `order` | 多项式阶数（1=线性, 2=二次） |
| `dim`   | 空间维度                     |

### 4.2 FiniteElementSpace — 串行有限元空间

```cpp
// 标量空间 (温度/电压)
auto scalar_fespace = std::make_unique<mfem::FiniteElementSpace>(mesh, fec);

// 矢量空间 (位移，vdim=dim)
auto vector_fespace = std::make_unique<mfem::FiniteElementSpace>(mesh, fec, dim);
```

| 接口                                          | 说明                            | 数据流角色      |
| --------------------------------------------- | ------------------------------- | --------------- |
| `GetVSize()`                                  | 返回 VSize                      | L-vector 的大小 |
| `GetTrueVSize()`                              | 串行时 = VSize                  | T-vector 的大小 |
| `GetEssentialTrueDofs(bdr_marker, tdof_list)` | 获取本质边界条件涉及的 DOF 列表 | 用于矩阵消元    |

### 4.3 ParFiniteElementSpace — 并行有限元空间

```cpp
auto par_fespace = std::make_unique<mfem::ParFiniteElementSpace>(pmesh, fec);
```

`ParFiniteElementSpace` 继承自 `FiniteElementSpace`，增加了并行 DOF 管理。

| 接口                        | 说明                                      | 数据流角色                 |
| --------------------------- | ----------------------------------------- | -------------------------- |
| `GetVSize()`                | 本地 L-vector 大小（含共享节点冗余）      | GridFunction 的 size       |
| `GetTrueVSize()`            | 本地 T-vector 大小（去冗余）              | 线性系统的本地维度         |
| `GlobalTrueVSize()`         | 全局总 DOF（所有 rank 的 TrueVSize 之和） | 仅用于日志                 |
| `GetProlongationMatrix()`   | 返回 P 矩阵：T-vector → L-vector          | **关键**：求解后恢复完整解 |
| `GetEssentialTrueDofs(...)` | T-vector 空间中的本质 DOF                 | 并行消元                   |

### 4.4 本项目的空间管理

所有有限元空间统一由 `FEConfig` 管理（[Config.hpp](include/fem/frontend/Config.hpp)），直接使用并行类型（parallel-only 架构）：

```cpp
struct FEConfig {
    std::unique_ptr<mfem::Mesh> serial_mesh;                    // 仅用于加载/构建
    std::unique_ptr<mfem::ParMesh> pmesh;                       // 并行网格
    std::unique_ptr<mfem::FiniteElementCollection> scalar_fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> scalar_fespace; // 标量并行空间
    std::unique_ptr<mfem::FiniteElementCollection> vector_fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> vector_fespace; // 矢量并行空间

    mfem::ParMesh &GetMesh() { return *pmesh; }
    mfem::ParFiniteElementSpace &GetScalarFESpace() { return *scalar_fespace; }
    mfem::ParFiniteElementSpace &GetVectorFESpace() { return *vector_fespace; }
};
```

后续所有 Solver 和 Assembler 直接使用 `ParFiniteElementSpace&`、`ParGridFunction`、`ParBilinearForm`、`ParLinearForm`，
无需 `dynamic_cast` 或串/并行分支判断。即使 `np=1`，同样走 `ParMesh` + `ParFiniteElementSpace` 路径。

---

## 5. GridFunction — 场变量载体

`GridFunction` 是 MFEM 中最核心的数据容器，它是一个绑定了有限元空间的 `Vector`。

### 5.1 本项目中的 GridFunction 实例

| 变量名            | 空间     | 含义            | 文件                                                                     |
| ----------------- | -------- | --------------- | ------------------------------------------------------------------------ |
| `temperature_`    | 标量空间 | 温度场          | [ThermalFieldSolver.cpp](src/physics/ThermalFieldSolver.cpp)             |
| `temperature_cn_` | 标量空间 | CN 方案的温度解 | [ThermalFieldSolver.cpp](src/physics/ThermalFieldSolver.cpp)             |
| `voltage_`        | 标量空间 | 电压场          | [ElectrostaticFieldSolver.cpp](src/physics/ElectrostaticFieldSolver.cpp) |
| `displacement_`   | 矢量空间 | 位移场          | [MechanicalFieldSolver.cpp](src/physics/MechanicalFieldSolver.cpp)       |
| `T_old`           | 标量空间 | 上一时间步温度  | [TransientSolver.cpp](src/transient/TransientSolver.cpp)                 |

### 5.2 关键接口

| 接口                                   | 说明                                      | 向量空间 |
| -------------------------------------- | ----------------------------------------- | -------- |
| `Size()` / 隐含大小                    | = VSize (L-vector)                        | L        |
| `GetTrueDofs(tv)`                      | 从 L-vector 提取 T-vector（restriction）  | L → T    |
| `SetFromTrueDofs(tv)`                  | 从 T-vector 恢复 L-vector（prolongation） | T → L    |
| `ProjectBdrCoefficient(coeff, marker)` | 将系数投影到边界节点                      | L        |
| `GetValue(T, ip)`                      | 在物理点处求值（标量）                    | L        |
| `GetVectorValue(T, ip, vec)`           | 在物理点处求值（矢量）                    | L        |
| `GetGradient(T, grad)`                 | 在物理点处求梯度                          | L        |
| `Min()` / `Max()`                      | 全局最小/最大值                           | L        |
| `operator=` (double)                   | 设置所有节点为常数                        | L        |

### 5.3 GridFunction 的内存布局（L-vector）

```
                        GridFunction (size = VSize)
    ┌──────────────────────────────────────────────────┐
    │  DOF_0  DOF_1  DOF_2  ... DOF_(VSize-1)          │
    │  ↑       ↑      ↑                                │
    │  本地    本地   共享节点的                         │
    │  内部    内部   冗余副本                           │
    │  节点    节点                                     │
    └──────────────────────────────────────────────────┘
    
    串行: 所有 DOF 都是"真实的"，没有冗余
    并行: 一些 DOF 是共享节点的冗余副本，它们的值
          通过 Prolongation/Restriction 与其他 rank 同步
```

---

## 6. 双线性形式与线性形式

### 6.1 创建方式 — 直接使用并行形式

本项目统一使用并行形式，各 Assembler 直接构造 `ParBilinearForm` / `ParLinearForm`：

```cpp
// 静电场组装 — ElectrostaticAssembler.cpp
auto bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
auto lf = std::make_unique<mfem::ParLinearForm>(&fespace);

// 热场组装 — ThermalAssembler.cpp
auto bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
auto lf = std::make_unique<mfem::ParLinearForm>(&fespace);

// 力学组装 — MechanicalAssembler.cpp
auto bf = std::make_unique<mfem::ParBilinearForm>(&fespace);
auto lf = std::make_unique<mfem::ParLinearForm>(&fespace);
```

相关文件：
- [ElectrostaticAssembler.hpp](include/fem/assembly/ElectrostaticAssembler.hpp) / [ElectrostaticAssembler.cpp](src/assembly/ElectrostaticAssembler.cpp)
- [ThermalAssembler.hpp](include/fem/assembly/ThermalAssembler.hpp) / [ThermalAssembler.cpp](src/assembly/ThermalAssembler.cpp)
- [MechanicalAssembler.hpp](include/fem/assembly/MechanicalAssembler.hpp) / [MechanicalAssembler.cpp](src/assembly/MechanicalAssembler.cpp)

### 6.2 BilinearForm / ParBilinearForm

**物理含义**：$ a(u, v) = \int_\Omega (\text{积分子}) \, d\Omega $ 对应的矩阵 $ K_{ij} $。

| 接口                                     | 说明                                         |
| ---------------------------------------- | -------------------------------------------- |
| `AddDomainIntegrator(integ)`             | 添加域积分子                                 |
| `AddBoundaryIntegrator(integ, marker)`   | 添加边界积分子（仅在 marker=1 的边界上积分） |
| `Assemble()`                             | 执行单元级组装 → 内部稀疏矩阵                |
| `Finalize()`                             | 完成稀疏矩阵的 CSR 格式化                    |
| `SpMat()`                                | 返回底层 `SparseMatrix&`（串行矩阵）         |
| `FormLinearSystem(tdofs, x, b, A, X, B)` | **核心**：消去本质 BC 后得到可求解的线性系统 |
| `RecoverFEMSolution(X, b, x)`            | 从求解结果恢复完整 FE 解                     |

**ParBilinearForm 独有**：

| 接口                 | 说明                                                               |
| -------------------- | ------------------------------------------------------------------ |
| `ParallelAssemble()` | 将串行 `SparseMatrix` 转为分布式 `HypreParMatrix`（涉及 MPI 通讯） |

### 6.3 LinearForm / ParLinearForm

**物理含义**：$ L(v) = \int_\Omega f \cdot v \, d\Omega $ 对应的右端向量 $ F_i $。

| 接口                                   | 说明                             |
| -------------------------------------- | -------------------------------- |
| `AddDomainIntegrator(integ)`           | 添加域积分子（源项）             |
| `AddBoundaryIntegrator(integ, marker)` | 添加边界积分子（通量/Robin RHS） |
| `Assemble()`                           | 执行单元级组装                   |
| `Size()`                               | = VSize (L-vector)               |

**ParLinearForm 独有**：

| 接口                  | 说明                                                |
| --------------------- | --------------------------------------------------- |
| `ParallelAssemble(F)` | 将 L-vector 的 RHS 收集(restrict)为 T-vector 的 RHS |

### 6.4 FormLinearSystem 详解

这是 MFEM 中最重要的接口之一。它同时处理：
1. 从双线性形式提取系统矩阵
2. 施加本质（Dirichlet）边界条件（行列消元）
3. 在并行模式下自动调用 `ParallelAssemble()`

```cpp
bf->FormLinearSystem(essential_tdofs, x_gf, rhs_lf, A, X, B);
// essential_tdofs: Array<int>        — T-vector 空间中的本质 DOF 列表
// x_gf:           ParGridFunction    — 初始猜测（含 BC 值），L-vector
// rhs_lf:         ParLinearForm      — 右端项，L-vector
// A:              OperatorPtr        — 输出：HypreParMatrix* (T-vector × T-vector)
// X, B:           Vector             — 输出：解和RHS，均为 TrueVSize 大小的 T-vector
```

**本项目使用点**：
- [ElectrostaticAssembler.cpp](src/assembly/ElectrostaticAssembler.cpp)：静电场求解
- [MechanicalFieldSolver.cpp](src/physics/MechanicalFieldSolver.cpp)：力学求解

> **注意**：热场求解器 **不使用** `FormLinearSystem`——
> 它手动构建 $ A = \alpha_k K + \alpha_c C $ 并手动消元，以便在时间步进中缓存矩阵。

---

## 7. 矩阵与向量

### 7.1 SparseMatrix — 本地稀疏矩阵

| 接口               | 说明                          | 本项目使用                           |
| ------------------ | ----------------------------- | ------------------------------------ |
| `Add(α, A, β, B)`  | 返回 $ \alpha A + \beta B $   | 构建 $ A = \alpha_k K + \alpha_c C $ |
| `AddMult(x, y, α)` | $ y \mathrel{+}= \alpha A x $ | RHS 构建                             |

### 7.2 HypreParMatrix — 并行分布式矩阵

MFEM 对 Hypre 的 `parcsr_matrix` 的封装。每个 rank 持有矩阵的若干行（行数 = TrueVSize）。

| 接口                       | 说明                              | 本项目使用                  |
| -------------------------- | --------------------------------- | --------------------------- |
| `Add(α, A, β, B)`          | 返回分布式 $ \alpha A + \beta B $ | 并行 thermal 的 $ A_{par} $ |
| `Mult(α, x, β, y)`         | $ y = \alpha Ax + \beta y $       | 并行 RHS 构建               |
| `EliminateRowsCols(tdofs)` | 并行消元                          | 并行 thermal 求解           |

### 7.3 Vector — 密集向量

| 接口                        | 说明                           |
| --------------------------- | ------------------------------ |
| `SetSize(n)`                | 分配 n 个元素                  |
| `operator=` (double)        | 设置所有元素为常数             |
| `Add(α, v)`                 | $ this \mathrel{+}= \alpha v $ |
| `operator*` (inner product) | $ \sum x_i y_i $               |
| `Size()`                    | 元素个数                       |
| `operator()(i)`             | 第 i 个元素（读写）            |

**本项目中向量的大小约定**：

| 场景                                 | 向量大小  | 空间     |
| ------------------------------------ | --------- | -------- |
| GridFunction                         | VSize     | L-vector |
| 线性系统 RHS `b`                     | TrueVSize | T-vector |
| 线性系统解 `x`                       | TrueVSize | T-vector |
| `FormLinearSystem` 输出 B, X         | TrueVSize | T-vector |
| ParLinearForm::ParallelAssemble 输出 | TrueVSize | T-vector |

### 7.4 Array\<int\>

主要用于边界标记数组：

```cpp
mfem::Array<int> marker(num_bdr);  // 长度 = 边界属性个数
marker = 0;                        // 全部标记为"不施加"
marker[attr - 1] = 1;             // 属性 attr 上施加 BC
```

### 7.5 OperatorPtr

`OperatorPtr` 是一个智能指针包装，自动管理 `FormLinearSystem` 返回的算子内存。它内部持有指向 `SparseMatrix` 或 `HypreParMatrix` 的指针。

---

## 8. 积分子（Integrators）

积分子是有限元组装的基本构件。每个积分子实现一个弱形式项的单元级矩阵/向量组装。

### 8.1 双线性形式积分子

| 积分子                                   | 弱形式                                                                                                                                        | 物理含义             | 使用场景                 |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------- | ------------------------ |
| `CustomDiffusionIntegrator(k)`           | $ \int_\Omega k \nabla u \cdot \nabla v \, d\Omega $                                                                                          | 热传导/电传导        | 热场 K 矩阵, 电场 K 矩阵 |
| `MassIntegrator(ρcₚ)`                    | $ \int_\Omega \rho c_p \, u \, v \, d\Omega $                                                                                                 | 热容（瞬态质量矩阵） | 热场 C 矩阵              |
| `ElasticityIntegrator(λ, μ)`             | $ \int_\Omega [\lambda (\nabla\cdot\mathbf{u})(\nabla\cdot\mathbf{v}) + 2\mu \, \varepsilon(\mathbf{u}):\varepsilon(\mathbf{v})] \, d\Omega $ | 线弹性刚度           | 力学 K 矩阵              |
| `CustomBoundaryMassIntegrator(l)`        | $ \int_{\Gamma_R} l \, u \, v \, d\Gamma $                                                                                                    | Robin BC 左端项      | 热场/电场 Robin BC       |
| `NormalDisplacementPenaltyIntegrator(β)` | $ \int_{\Gamma} \beta \, (u_n)(v_n) \, d\Gamma $                                                                                              | 法向位移惩罚         | 力学法向约束             |

### 8.2 线性形式积分子

| 积分子                              | 弱形式                                                                   | 物理含义        | 使用场景      |
| ----------------------------------- | ------------------------------------------------------------------------ | --------------- | ------------- |
| `CustomDomainLFIntegrator(f)`       | $ \int_\Omega f \, v \, d\Omega $                                        | 体源项          | 热源/电流源   |
| `VectorDomainLFIntegrator(f)`       | $ \int_\Omega \mathbf{f} \cdot \mathbf{v} \, d\Omega $                   | 体力            | 力学体力      |
| `ThermalStrainLFIntegrator(...)`    | $ \int_\Omega \sigma^{th}_{ij} \varepsilon_{ij}(\mathbf{v}) \, d\Omega $ | 热应变等效载荷  | 热-力耦合     |
| `CustomBoundaryLFIntegrator(q)`     | $ \int_{\Gamma_R} q \, v \, d\Gamma $                                    | Robin BC 右端项 | 热场 Robin BC |
| `VectorBoundaryFluxLFIntegrator(p)` | $ \int_\Gamma p \, \mathbf{n} \cdot \mathbf{v} \, d\Gamma $              | 压力/面力       | 力学压力 BC   |

### 8.3 自定义积分子结构

所有自定义积分子继承自 `BilinearFormIntegrator` 或 `LinearFormIntegrator`，需要实现：

```cpp
// 双线性形式：计算单元矩阵 K_e
void AssembleElementMatrix(const FiniteElement &el, 
                           ElementTransformation &trans,
                           DenseMatrix &elmat) override;

// 线性形式：计算单元向量 f_e
void AssembleRHSElementVect(const FiniteElement &el,
                            ElementTransformation &trans,
                            Vector &elvect) override;
```

组装过程中 MFEM 遍历所有单元，对每个单元调用上述虚函数，
将得到的 `elmat` / `elvect` 散列到全局矩阵/向量中。

---

## 9. 系数（Coefficients）

系数表示 PDE 中的物理参数（导热率、电导率等），在积分点处求值。

### 9.1 MFEM 内建系数

| 类                               | 功能                 | 本项目用途                |
| -------------------------------- | -------------------- | ------------------------- |
| `ConstantCoefficient(val)`       | 常数                 | Dirichlet BC 值、默认源项 |
| `PWConstCoefficient(vals)`       | 按域属性分段常数     | 多域源项                  |
| `GridFunctionCoefficient(gf)`    | 从 GridFunction 插值 | 温度 → 力学耦合           |
| `VectorConstantCoefficient(vec)` | 常矢量               | 均匀体力                  |
| `VectorArrayCoefficient(dim)`    | 按分量组合矢量系数   | 多域体力                  |

### 9.2 本项目自定义系数

| 类                             | 功能                                          | 文件                                                       |
| ------------------------------ | --------------------------------------------- | ---------------------------------------------------------- |
| `PiecewiseConstantCoefficient` | 根据材料数据库按域赋值，预评估常数            | [CoefficientManager.cpp](src/coeff/CoefficientManager.cpp) |
| `ExpressionCoefficient`        | 通过 muparser 支持温度依赖表达式 `"sigma(T)"` | [CoefficientManager.cpp](src/coeff/CoefficientManager.cpp) |
| `JouleHeatingCoefficient`      | 计算 $ \sigma                                 | \nabla V                                                   | ^2 $ | [CoefficientManager.cpp](src/coeff/CoefficientManager.cpp) |

所有系数的核心虚函数：
```cpp
double Eval(ElementTransformation &T, const IntegrationPoint &ip) override;
```
在积分点 `(T, ip)` 处返回系数值。

---

## 10. 线性求解器

### 10.1 求解器接口

```cpp
class ILinearSolver {
    virtual void Solve(Operator &A, const Vector &b, Vector &x) = 0;  // 首次求解（含分解）
    virtual void Mult(const Vector &b, Vector &x) = 0;                // 复用分解
};
```

### 10.2 可用求解器

| 类                  | 底层                     | 适用模式  | 特性                |
| ------------------- | ------------------------ | --------- | ------------------- |
| `MfemPcgSolver`     | MFEM CGSolver + 预条件器 | 串行/并行 | 迭代，对称正定      |
| `MfemPardisoSolver` | MKL PARDISO              | 串行      | 直接法，LU 分解缓存 |
| `MfemMumpsSolver`   | MUMPS                    | 串行/并行 | 直接法，分布式分解  |
| `MfemUmfpackSolver` | UMFPACK                  | 串行      | 直接法              |

### 10.3 求解器工厂

```cpp
auto solver = solver::CreateLinearSolver("pardiso", rel_tol, abs_tol, max_iter, print_level);
```

### 10.4 缓存策略

本项目的求解器在首次调用 `Solve()` 时进行矩阵分解，后续时间步只需 `Mult()` 即可复用分解，
直到矩阵变化（如 dt 改变导致 $ A = \alpha K + \beta C $ 重构）才重新 `Solve()`。

---

## 11. 边界条件处理

### 11.1 本质（Dirichlet）边界条件

**流程**：
1. 构建 marker 数组：`Array<int> marker(num_bdr); marker[attr-1] = 1;`
2. 获取受影响的 DOF：`fespace.GetEssentialTrueDofs(marker, tdof_list)`
3. 将 BC 值投影到解向量：`gf.ProjectBdrCoefficient(coeff, marker)`
4. 消元：
   - 使用 `FormLinearSystem`（自动处理）—— 静电/力学
   - 手动构建 $ A $ 后用 `EliminateRowsCols`—— 热场

### 11.2 Robin 边界条件

不修改 DOF 列表，而是通过积分子添加弱形式项：

$$ \text{LHS: } \int_{\Gamma_R} l \cdot u \cdot v \, d\Gamma, \quad \text{RHS: } \int_{\Gamma_R} q \cdot v \, d\Gamma $$

### 11.3 惩罚法边界条件（力学法向位移）

通过大惩罚系数 $ \beta \sim 10^{12} $ 近似约束法向位移为零：

$$ \int_\Gamma \beta \, (u_n)(v_n) \, d\Gamma \approx 0 \Rightarrow u_n \approx 0 $$

---

## 12. 完整求解数据流

### 12.1 热场稳态求解（手动矩阵路径）

```
[组装阶段]
                          ┌──────────────────────────┐
  ParBilinearForm         │ AddDomainIntegrator(k)   │
                          │ AddBoundaryIntegrator(l)  │
                          │ Assemble() → SparseMatrix │
                          │ Finalize()                │
                          │ SpMat() → K_sparse        │
                          │ ParallelAssemble()→K_par  │
                          └──────────────────────────┘
                          
  ParLinearForm           ┌──────────────────────────┐
                          │ AddDomainIntegrator(f)    │
                          │ AddBoundaryIntegrator(q)  │
                          │ Assemble()                │
                          │ ParallelAssemble(F)       │
                          └──────────────────────────┘

[求解阶段]

  K_par (TrueVSize × TrueVSize per rank, HypreParMatrix)
  EliminateRowsCols(tdofs)
  → A_par
  
  b = F (TrueVSize)                           ← ParallelAssemble 后已是 T-vector
  temperature_.GetTrueDofs(x_true)             ← L → T
  ApplyBCToRHS: b -= K_par * bc_vec; b[td] = bc_val
  
  solver->Solve(A_par, b, x_true)             // x_true: TrueVSize
  
  P->Mult(x_true, temperature_)               ← T → L (prolongation)
  // temperature_: VSize，可用于后处理
```

### 12.2 热场瞬态求解（BE+CN）

```
BE (Backward Euler):
  A = (1/dt)C + K
  b = (1/dt)C·T_old + F

CN (Crank-Nicolson):
  A = (1/dt)C + 0.5K
  b = (1/dt)C·T_old - 0.5K·T_old + 0.5(F_prev + F_current)

向量空间（parallel-only）:
  T_old_true ← T_old_.GetTrueDofs()      // L → T
  C_par->Mult(inv_dt, T_old_true, 0.0, b)  // 在 T 空间计算
  solver->Solve(A_par, b, x_true)
  P->Mult(x_true, temperature_)           // T → L
```

### 12.3 静电场求解（FormLinearSystem 路径）

```
  bf->Assemble()   // ParBilinearForm
  lf->Assemble()   // ParLinearForm
  
  fespace.GetEssentialTrueDofs(bdr, tdofs)
  voltage_.ProjectBdrCoefficient(...)
  
  bf->FormLinearSystem(tdofs, voltage_, *lf, A, X, B)
  // FormLinearSystem 内部：
  //   ParallelAssemble → HypreParMatrix, 消元 → A (HypreParMatrix*)
  //   X, B 均为 TrueVSize 大小的 T-vector
  
  solver->Solve(*A, B, X)
  
  bf->RecoverFEMSolution(X, *lf, voltage_)
  // RecoverFEMSolution 内部：
  //   反向消元 + P->Mult → 恢复完整 L-vector (VSize)
```

### 12.4 力学场求解（缓存刚度 + RHS-only）

```
Phase 1: CacheStiffnessMatrix()
  AssembleStiffness(fespace, λ, μ, penalty_coeffs, nd_markers)
  → cached_bilinear_ (ParBilinearForm)

Phase 2: SolveRHSOnly() (每个时间步调用)
  lf = AssembleRHS(fespace, body_force, pressure, thermal_strain, ...)
  displacement_.ProjectBdrCoefficient(...)
  
  cached_bilinear_->FormLinearSystem(tdofs, displacement_, *lf, A, X, B)
  
  if (!factorized)
      solver->Solve(*A, B, X)   // 首次分解
  else
      solver->Mult(B, X)        // 复用分解
  
  cached_bilinear_->RecoverFEMSolution(X, *lf, displacement_)
```

### 12.5 瞬态自适应误差估计中的 MPI 通讯

```
WRMS 误差计算 (ComputeAdaptiveError):

  基于 T-vector（TrueVSize）计算，避免共享节点冗余影响:
    local_sum_sq = Σ [(T_cn_i - T_be_i) / (abstol + reltol·|T_cn_i|)]²
    local_N = T_be.Size()   // = TrueVSize

  mfem::InnerProduct(comm, ...): 全局汇总
    global_sum_sq = Σ(all ranks) local_sum_sq
    global_N = Σ(all ranks) local_N

  WRMS = √(global_sum_sq / global_N)

Picard 收敛判据 (PicardCoupler::RunPicard):

  基于 T-vector（TrueVSize）计算:
    T_new.GetTrueDofs(true_vec)
    diff_sq = mfem::InnerProduct(comm, delta, delta)
    ref_sq  = mfem::InnerProduct(comm, T_prev_true, T_prev_true)
  
  rel_change = √(diff_sq) / √(ref_sq)
```

---

## 13. 输出与导出

### 13.1 ParaViewDataCollection

通过 [ParaviewExporter.cpp](src/output/ParaviewExporter.cpp) 导出 VTK 文件，支持并行：

```cpp
mfem::ParaViewDataCollection collection(name, &mesh);
collection.SetPrefixPath(output_dir);
collection.SetLevelsOfDetail(1);
collection.SetCycle(cycle);           // 时间步索引
collection.SetTime(time);            // 物理时间
collection.RegisterField("T", &gf);  // 注册 GridFunction
collection.SetHighOrderOutput(true);
collection.SetDataFormat(mfem::VTKFormat::BINARY);
collection.Save();
```

并行时 `ParaViewDataCollection` 自动识别 `ParMesh`，各 rank 写出自己分区的 `.vtu` 文件，
rank 0 额外写出 `.pvtu` 汇总文件。ParaView 可直接打开 `.pvtu` 查看完整场。

### 13.2 文本导出

[SolutionTextExporter.cpp](src/output/SolutionTextExporter.cpp) 接收串行化后的 `Mesh&` 和 `GridFunction&`，通过逐顶点求值输出节点数据。

在并行模式下，[ResultExporter.cpp](src/output/ResultExporter.cpp) 先通过
`ParMesh::GetSerialMesh(0)` + `ParGridFunction::GetSerialGridFunction(0, serial_mesh)` 
将并行数据收集到 rank 0 的串行对象上，再调用 `SolutionTextExporter` 导出。
这使得 txt 导出在任意 `np` 下均可工作。

关键 MFEM 接口调用链：
```
// ResultExporter 中的并行→串行收集
auto serial_mesh = mesh_.GetSerialMesh(0);
auto serial_gf = gf.GetSerialGridFunction(0, serial_mesh);

// SolutionTextExporter 中的串行网格遍历
mesh.GetElementVertices(e, verts)          → 获取顶点 → 单元映射
mesh.GetElementTransformation(e)           → 参考 → 物理坐标变换
Geometries.GetVertices(geom).IntPoint(lv)  → 顶点在参考空间的积分点
scalar_field.GetValue(T, ip)               → 在物理点插值求值
vector_field.GetVectorValue(T, ip, vec)    → 矢量场插值
```

### 13.3 并行输出中的 rank 协调

```cpp
// ResultExporter 构造函数中（直接使用 ParMesh）:
int rank, np;
MPI_Comm_rank(mesh.GetComm(), &rank);
MPI_Comm_size(mesh.GetComm(), &np);

if (rank == 0)
    std::filesystem::remove_all(output_dir_);  // 仅 rank 0 清理
if (np > 1)
    MPI_Barrier(mesh.GetComm());               // 等待清理完成
std::filesystem::create_directories(output_dir_);   // 所有 rank 创建目录
```

---

## 14. MPI 通讯点汇总

下表列出本项目中所有发生 MPI 通讯的位置及其通讯方式：

| 位置                                  | 通讯类型      | 目的                                  | 数据                          |
| ------------------------------------- | ------------- | ------------------------------------- | ----------------------------- |
| `ParMesh(comm, mesh)`                 | 集合通讯      | METIS 分区 + 分发子网格               | 整个网格                      |
| `ParFiniteElementSpace` 构造          | 集合通讯      | 协调共享 DOF 编号                     | DOF → 全局映射                |
| `ParBilinearForm::ParallelAssemble()` | 点对点 + 集合 | 共享 DOF 的矩阵行合并                 | SparseMatrix → HypreParMatrix |
| `ParLinearForm::ParallelAssemble(F)`  | 点对点        | 共享 DOF 的 RHS 值合并                | L-vector → T-vector           |
| `HypreParMatrix::Mult()`              | 点对点        | 并行矩阵-向量乘法（需要邻居的解值）   | T-vector                      |
| `HypreParMatrix::Add()`               | 本地          | 矩阵加法（各 rank 独立，无通讯）      | —                             |
| `HypreParMatrix::EliminateRowsCols()` | 本地          | BC 消元（各 rank 独立）               | —                             |
| `GetProlongationMatrix()->Mult()`     | 点对点        | T → L 展开（需接收邻居的共享 DOF 值） | T-vector → L-vector           |
| `GetTrueDofs()`                       | 本地          | L → T 抽取（纯本地操作）              | —                             |
| `mfem::InnerProduct` (WRMS)           | 集合通讯      | 全局误差汇总                          | TrueVSize vector              |
| `mfem::InnerProduct` (Picard)         | 集合通讯      | 全局收敛判据                          | TrueVSize vector              |
| `GetSerialMesh/GetSerialGridFunction` | 集合通讯      | 并行数据收集到 rank 0 用于 txt 输出   | 全局网格 + 场数据             |
| `MPI_Barrier` (输出)                  | 集合通讯      | 同步文件系统操作                      | —                             |
| `ParaViewDataCollection::Save()`      | 文件 I/O      | 各 rank 写各自 .vtu，rank 0 写 .pvtu  | 场数据                        |
| Hypre/MUMPS 求解器内部                | 复杂通讯      | 分布式 LU 分解 / Krylov 迭代          | 矩阵/向量分块                 |

### 通讯频率分析

- **一次性（初始化）**：ParMesh、ParFESpace、ParallelAssemble（K, C 矩阵）
- **每时间步**：ParallelAssemble（RHS F）、矩阵-向量乘法、求解器、Prolongation
- **仅在自适应/Picard 中**：MPI_Allreduce (少量 double)

这意味着并行开销主要来自：
1. 求解器内部通讯（Hypre/MUMPS 的分布式分解/迭代）
2. 每步的 RHS ParallelAssemble + 矩阵-向量乘 + Prolongation

对于本项目的中等规模网格（数千 DOF），并行加速比受限于通讯-计算比，
通常在 4-8 核时达到最佳效率。更大的网格才能充分发挥并行优势。

---

## 附录 A：本项目并行代码路径示例

> 本项目统一使用 parallel-only 架构（即使 `np=1` 也使用 `ParMesh` + `ParFiniteElementSpace`），
> 不存在串行分支。以下展示本项目实际使用的并行代码路径。

```cpp
// ============= 组装 =============
ParBilinearForm bf(par_fespace);
bf.Assemble();
HypreParMatrix *K_par = bf.ParallelAssemble();  // MPI 通讯点

// ============= RHS =============
ParLinearForm lf(par_fespace);
lf.Assemble();  // L-vector
Vector F(par_fespace->GetTrueVSize());
lf.ParallelAssemble(F);  // L → T, MPI 通讯点

// ============= 求解 (手动矩阵路径，如热场) =============
K_par->EliminateRowsCols(tdofs);
gf.GetTrueDofs(x_true);  // L → T (本地)
solver.Solve(*K_par, F, x_true);  // x_true: TrueVSize, MPI 通讯
P->Mult(x_true, gf);  // T → L, MPI 通讯点

// ============= 求解 (FormLinearSystem 路径，如静电/力学) =============
bf.FormLinearSystem(tdofs, gf, lf, A, X, B);  // 自动 ParallelAssemble + 消元
solver.Solve(*A, B, X);
bf.RecoverFEMSolution(X, lf, gf);             // 自动 T → L

// ============= 结果导出（并行 → 串行收集） =============
auto serial_mesh = pmesh.GetSerialMesh(0);
auto serial_gf = pgf.GetSerialGridFunction(0, serial_mesh);
if (rank == 0)
    SolutionTextExporter::ExportScalarNodalTxt(..., serial_mesh, serial_gf, ...);
```

## 附录 B：本项目 MFEM 接口完整索引

以下按字母顺序列出所有使用到的 MFEM 类和函数：

### 类

| 类名                             | 头文件     | 用途                                  |
| -------------------------------- | ---------- | ------------------------------------- |
| `Array<int>`                     | `mfem.hpp` | 边界标记、DOF 列表                    |
| `BilinearFormIntegrator`         | `mfem.hpp` | 自定义积分子基类                      |
| `CGSolver`                       | `mfem.hpp` | CG 迭代求解器                         |
| `ConstantCoefficient`            | `mfem.hpp` | 常数系数                              |
| `DenseMatrix`                    | `mfem.hpp` | 单元矩阵                              |
| `DiffusionIntegrator`            | `mfem.hpp` | 扩散积分子（被自定义版替代）          |
| `ElasticityIntegrator`           | `mfem.hpp` | 弹性积分子                            |
| `ElementTransformation`          | `mfem.hpp` | 坐标变换                              |
| `FiniteElementSpace`             | `mfem.hpp` | FE 空间基类（用于串行化后的数据导出） |
| `GridFunction`                   | `mfem.hpp` | 场变量基类                            |
| `GridFunctionCoefficient`        | `mfem.hpp` | GF → 系数                             |
| `H1_FECollection`                | `mfem.hpp` | H1 有限元族                           |
| `Hypre`                          | `mfem.hpp` | Hypre 初始化                          |
| `HypreParMatrix`                 | `mfem.hpp` | 并行稀疏矩阵                          |
| `IntegrationPoint`               | `mfem.hpp` | 积分点                                |
| `IntegrationRule`                | `mfem.hpp` | 积分规则                              |
| `LinearFormIntegrator`           | `mfem.hpp` | 自定义积分子基类                      |
| `MassIntegrator`                 | `mfem.hpp` | 质量积分子                            |
| `Mesh`                           | `mfem.hpp` | 串行网格（用于加载和串行化导出）      |
| `Mpi`                            | `mfem.hpp` | MPI 工具类                            |
| `OperatorPtr`                    | `mfem.hpp` | 算子智能指针                          |
| `PWConstCoefficient`             | `mfem.hpp` | 分段常数系数                          |
| `ParBilinearForm`                | `mfem.hpp` | 并行双线性形式（本项目主要使用）      |
| `ParFiniteElementSpace`          | `mfem.hpp` | 并行 FE 空间（本项目主要使用）        |
| `ParGridFunction`                | `mfem.hpp` | 并行场变量（本项目主要使用）          |
| `ParLinearForm`                  | `mfem.hpp` | 并行线性形式（本项目主要使用）        |
| `ParMesh`                        | `mfem.hpp` | 并行网格（本项目主要使用）            |
| `ParaViewDataCollection`         | `mfem.hpp` | ParaView 输出                         |
| `Solver`                         | `mfem.hpp` | 求解器基类                            |
| `SparseMatrix`                   | `mfem.hpp` | 本地稀疏矩阵（ParallelAssemble 前）   |
| `VTKFormat`                      | `mfem.hpp` | VTK 格式枚举                          |
| `Vector`                         | `mfem.hpp` | 密集向量                              |
| `VectorArrayCoefficient`         | `mfem.hpp` | 分量组合矢量系数                      |
| `VectorBoundaryFluxLFIntegrator` | `mfem.hpp` | 边界通量积分子                        |
| `VectorConstantCoefficient`      | `mfem.hpp` | 常矢量系数                            |
| `VectorDomainLFIntegrator`       | `mfem.hpp` | 矢量域积分子                          |

### 全局函数 / 静态方法

| 函数                                       | 用途                               |
| ------------------------------------------ | ---------------------------------- |
| `mfem::Add(α, A, β, B)`                    | 矩阵线性组合                       |
| `mfem::add(α, x, β, y, z)`                 | 向量线性组合                       |
| `mfem::Geometries.GetVertices(geom)`       | 参考单元顶点                       |
| `mfem::Mpi::Init()`                        | MPI 初始化                         |
| `mfem::Mpi::IsInitialized()`               | 检查 MPI 状态                      |
| `mfem::Mpi::WorldSize()`                   | MPI 进程总数                       |
| `mfem::Hypre::Init()`                      | Hypre 初始化                       |
| `mfem::InnerProduct(comm, x, y)`           | MPI 全局内积（Picard/WRMS 中使用） |
| `MPI_Barrier()`                            | 同步屏障                           |
| `MPI_Comm_rank()`                          | 获取当前 rank 编号                 |
| `ParMesh::GetSerialMesh()`                 | 并行网格收集为串行网格（rank 0）   |
| `ParGridFunction::GetSerialGridFunction()` | 并行场收集为串行场（rank 0）       |

# MFEM MPI 接口与多进程架构

本文档整理 MFEM 并行求解常用接口和多进程数据流。普通 MFEM 对象见 [mfem_common_interfaces.md](mfem_common_interfaces.md)，本项目实际调用索引见 [mfem_api.md](mfem_api.md)，rank 内线程化分析见 [hybrid_parallel_mfem.md](hybrid_parallel_mfem.md)。

## 1. 并行程序骨架

典型入口：

```cpp
int main(int argc, char *argv[])
{
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();

    // build ParMesh, ParFiniteElementSpace, ParGridFunction...

    return 0;
}
```

本项目入口在 `main.cpp` 中初始化 MPI 和 Hypre，然后所有求解路径都使用并行对象；即使 `-np 1`，也运行同一套并行代码。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `mfem::Mpi::Init(argc, argv)` | 初始化 MPI。 |
| `mfem::Mpi::WorldRank()` | 当前 rank。 |
| `mfem::Mpi::WorldSize()` | 总 rank 数。 |
| `mfem::Hypre::Init()` | 初始化 Hypre 并行线性代数环境。 |
| `MPI_Comm_rank(comm, &rank)` | MPI 原生 rank 查询。 |
| `MPI_Comm_size(comm, &size)` | MPI 原生 size 查询。 |

## 2. 多进程对象关系

并行对象链路：

```text
Mesh on every rank
  -> ParMesh partitions ownership
  -> ParFiniteElementSpace defines local/true/global DOFs
  -> ParGridFunction stores rank-local L-vector
  -> ParBilinearForm / ParLinearForm assemble local contributions
  -> ParallelAssemble builds global distributed system
  -> HypreParMatrix + Vector true DOF solve
  -> Distribute or SetFromTrueDofs recovers field
```

| 串行概念 | 并行概念 | 说明 |
| -------- | -------- | ---- |
| `Mesh` | `ParMesh` | 分布式网格，每个 rank 拥有一部分单元和共享边界信息。 |
| `FiniteElementSpace` | `ParFiniteElementSpace` | 分布式有限元空间。 |
| `GridFunction` | `ParGridFunction` | 每个 rank 的本地场向量。 |
| `BilinearForm` | `ParBilinearForm` | 本地装配 + 并行组装接口。 |
| `LinearForm` | `ParLinearForm` | 本地 RHS 装配 + 并行组装接口。 |
| `SparseMatrix` | `HypreParMatrix` | 分布式稀疏矩阵。 |

## 3. ParMesh

### 3.1 构造

常见构造：

```cpp
mfem::Mesh serial_mesh("mesh.msh", 1, 1);
mfem::ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);
```

重要语义：

- 串行网格通常在每个 rank 上读入，然后由 `ParMesh` 重新划分。
- `ParMesh` 持有 rank 本地单元、共享面的通信信息和边界属性。
- 后续有限元空间必须基于 `ParMesh` 创建。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `Dimension()` | 空间维度。 |
| `GetNE()` | 当前 rank 的本地单元数。 |
| `GetNBE()` | 当前 rank 的本地边界单元数。 |
| `GetNV()` | 当前 rank 可见顶点数。 |
| `GetComm()` | 关联的 MPI communicator。 |
| `GetMyRank()` | 当前 rank。 |
| `GetNRanks()` | communicator 中 rank 数。 |
| `UniformRefinement()` | 并行均匀细化。 |
| `GetSerialMesh(rank)` | 收集/生成串行网格，常用于 rank 0 输出。 |

## 4. ParFiniteElementSpace

构造：

```cpp
mfem::H1_FECollection fec(order, dim);
mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
mfem::ParFiniteElementSpace vpfes(&pmesh, &fec, dim);
```

### 4.1 三类 DOF

| 概念 | 接口 | 说明 |
| ---- | ---- | ---- |
| Local DOF / L-vector | `GetVSize()` | rank 本地向量大小，包含共享自由度的本地副本。 |
| True DOF / T-vector | `GetTrueVSize()` | rank 拥有的真实自由度数量，用于线性系统。 |
| Global true DOF | `GlobalTrueVSize()` | 全局真实自由度总数。 |

直觉：

- `ParGridFunction` 使用 L-vector。
- `HypreParMatrix` 和线性求解器使用 true DOF 向量。
- 边界条件最终也要转换成 true DOF 编号。

### 4.2 常用接口

| 接口 | 说明 |
| ---- | ---- |
| `GetVSize()` | rank 本地 L-vector 长度。 |
| `GetTrueVSize()` | rank 本地 true DOF 长度。 |
| `GlobalTrueVSize()` | 全局 true DOF 总数。 |
| `GetEssentialTrueDofs(marker, tdofs)` | 根据边界 marker 得到本质边界 true DOF。 |
| `GetRestrictionMatrix()` | L-vector 到 true vector 的限制矩阵。 |
| `GetProlongationMatrix()` | true vector 到 L-vector 的延拓矩阵。 |
| `GetComm()` | communicator。 |
| `GetParMesh()` | 关联 `ParMesh`。 |

## 5. ParGridFunction

`ParGridFunction` 保存 rank 本地场变量。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `ParGridFunction(pfes)` | 在并行空间上创建场。 |
| `ProjectCoefficient(coeff)` | 投影系数到场。 |
| `GetTrueDofs(vec)` | 从 L-vector 抽取 true DOF 向量。 |
| `SetFromTrueDofs(vec)` | 用 true DOF 向量设置本地场。 |
| `Distribute(vec)` | 将 true DOF 解分发到 L-vector。 |
| `GetValue(elem, ip)` | 在 rank 本地单元上求值。 |
| `GetGradient(elem, ip, grad)` | 标量场梯度。 |
| `GetSerialGridFunction(rank)` | 收集成串行场，常用于 rank 0 文本输出。 |

使用要点：

- `GetValue()`、`GetGradient()` 只对本 rank 可见单元有意义。
- 线性求解后必须把 true DOF 解分发回 `ParGridFunction`，否则输出场不会更新。
- 收集到 rank 0 的串行场只适合输出和后处理，不应再参与并行求解。

## 6. 并行弱形式

### 6.1 ParBilinearForm

常见流程：

```cpp
mfem::ParBilinearForm a(&pfes);
a.AddDomainIntegrator(new mfem::DiffusionIntegrator(k));
a.Assemble();
a.Finalize();
std::unique_ptr<mfem::HypreParMatrix> A(a.ParallelAssemble());
```

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `AddDomainIntegrator()` | 添加体积分子。 |
| `AddBoundaryIntegrator()` | 添加边界积分子。 |
| `Assemble()` | 本地装配。 |
| `Finalize()` | 完成本地矩阵。 |
| `ParallelAssemble()` | 组装为分布式矩阵。 |
| `FormLinearSystem(ess_tdofs, x, b, A, X, B)` | 生成带本质边界处理的并行线性系统。 |
| `RecoverFEMSolution(X, b, x)` | 从 true DOF 解恢复有限元场。 |

### 6.2 ParLinearForm

常见流程：

```cpp
mfem::ParLinearForm b(&pfes);
b.AddDomainIntegrator(new mfem::DomainLFIntegrator(f));
b.Assemble();
mfem::Vector B;
b.ParallelAssemble(B);
```

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `AddDomainIntegrator()` | 添加体源。 |
| `AddBoundaryIntegrator()` | 添加边界源。 |
| `Assemble()` | 本地 RHS 装配。 |
| `ParallelAssemble(vec)` | 组装 true DOF RHS。 |

## 7. HypreParMatrix 与 true DOF 线性系统

`HypreParMatrix` 是 MFEM 中最常见的 MPI 分布式稀疏矩阵。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `Height()` / `Width()` | 本地算子尺寸。 |
| `GetGlobalNumRows()` | 全局行数。 |
| `Mult(x, y)` | 分布式矩阵向量乘。 |
| `EliminateRowsCols(ess_tdofs, diag)` | 消元本质边界条件。 |
| `Read_IJMatrix()` | 获取底层 IJ 矩阵指针，通常只在高级场景使用。 |

常见求解器组合：

```cpp
mfem::HypreBoomerAMG amg(*A);

mfem::CGSolver cg(comm);
cg.SetPreconditioner(amg);
cg.SetOperator(*A);
cg.Mult(B, X);
```

| 求解器 | 说明 |
| ------ | ---- |
| `CGSolver(comm)` | 并行 CG。 |
| `GMRESSolver(comm)` | 并行 GMRES。 |
| `HypreBoomerAMG` | AMG 预条件器。 |
| `HyprePCG` | Hypre PCG 包装。 |
| `MUMPSSolver` | 并行直接求解器，需要构建支持。 |
| `SuperLUSolver` | SuperLU_DIST 并行直接求解器，需要构建支持。 |
| `STRUMPACKSolver` | STRUMPACK 并行求解器，需要构建支持。 |
| `CPardisoSolver` | MKL Cluster PARDISO，需要构建支持。 |

## 8. 边界条件的并行处理

并行路径中，边界属性 marker 仍按属性编号构造：

```cpp
mfem::Array<int> marker(pmesh.bdr_attributes.Max());
marker = 0;
marker[attr - 1] = 1;
```

然后转换为 true DOF：

```cpp
mfem::Array<int> ess_tdofs;
pfes.GetEssentialTrueDofs(marker, ess_tdofs);
```

两种常见处理方式：

| 方式 | 说明 |
| ---- | ---- |
| `FormLinearSystem()` | 让 MFEM 统一处理本质边界、矩阵/RHS 和 true DOF 变量。 |
| 手动消元 | 先 `ParallelAssemble()`，再对 `HypreParMatrix` 与 RHS 做消元。适合需要缓存矩阵的流程。 |

本项目两种方式都有：静电场偏向 `FormLinearSystem()` 风格，热/力部分为了缓存矩阵会手动控制更多步骤。

## 9. 并行数据收集与输出

常见输出方案：

| 方案 | 接口 | 说明 |
| ---- | ---- | ---- |
| 分布式 ParaView | `ParaViewDataCollection` | 每个 rank 写自己的 piece，适合可视化。 |
| rank 0 串行文本 | `GetSerialMesh()` + `GetSerialGridFunction()` | 收集到 rank 0 后写 COMSOL 风格文本。 |
| 自定义 MPI gather | MPI 原生 gather/allgather | 适合专门的统计或抽样输出。 |

本项目文本导出流程：

```text
每个 rank 持有 ParGridFunction
  -> 收集为 rank 0 上的串行 Mesh/GridFunction
  -> rank 0 遍历节点并写 txt
```

注意：

- 分布式场的本地节点顺序不等于全局文本输出顺序。
- 做坐标匹配比较时，应以 `x,y,z` 为键，而不是依赖行号。
- rank 0 文本导出是后处理步骤，不参与求解通信。

## 10. 通信发生在哪里

常见通信点：

| 场景 | 典型接口 | 通信含义 |
| ---- | -------- | -------- |
| 构造 `ParMesh` | `ParMesh(comm, mesh)` | 网格分区和共享实体建立。 |
| 构造并行空间 | `ParFiniteElementSpace` | 建立共享 DOF 与 true DOF 映射。 |
| 并行装配 | `ParallelAssemble()` | 汇总共享自由度贡献。 |
| 形成线性系统 | `FormLinearSystem()` | 本质边界和 true DOF 系统构造。 |
| Krylov 求解 | `CGSolver::Mult()` 等 | 矩阵向量乘、点积、预条件通信。 |
| 场恢复 | `Distribute()` / `RecoverFEMSolution()` | true DOF 到本地 L-vector 同步。 |
| 全局范数 | `InnerProduct(comm, x, y)` 等 | reduce 操作。 |
| 文本输出 | `GetSerialGridFunction()` | 收集到指定 rank。 |

## 11. 多进程架构建议

### 11.1 推荐主路径

```text
每个 rank:
  读配置
  读串行网格
  构造 ParMesh
  构造 ParFiniteElementSpace
  本地装配
  并行组装 true DOF 系统
  并行求解
  本地恢复 ParGridFunction

rank 0:
  额外执行文本输出和 Python 对比脚本
```

### 11.2 避免的常见错误

| 错误 | 后果 |
| ---- | ---- |
| 用 `GetVSize()` 创建 true DOF RHS | 维度错配或隐式错误。 |
| 求解后忘记 `Distribute()` / `RecoverFEMSolution()` | 输出场仍是旧值。 |
| 多个 rank 同时写同一个文本文件 | 文件竞争或内容交错。 |
| 用行号比较分布式输出 | 并行分区变化后顺序不稳定。 |
| 在 rank 0 串行对象上继续并行求解 | 破坏并行数据一致性。 |

## 12. 与本项目当前实现的对应关系

| 项目模块 | 并行职责 |
| -------- | -------- |
| `ConfigLoader::BuildFE()` | 读取网格，构造 `ParMesh`、标量/矢量并行空间。 |
| `ElectrostaticFieldSolver` | 装配静电并行系统，求解电势 true DOF。 |
| `ThermalFieldSolver` | 缓存热传导刚度/质量矩阵，按时间步求解温度 true DOF。 |
| `MechanicalFieldSolver` | 缓存力学刚度矩阵，按快照求解位移 true DOF。 |
| `PicardCoupler` | 在 true DOF 上计算温度变化范数并判断收敛。 |
| `ResultExporter` | 分布式 ParaView 输出和 rank 0 文本导出。 |
| `Application` | 求解后只在 rank 0 调用 Python 对比脚本。 |

当前项目已经把“MPI 多进程”作为默认求解骨架；下一步若要做 rank 内线程化，应优先保持这些 MPI 对象的所有权和通信语义不变。

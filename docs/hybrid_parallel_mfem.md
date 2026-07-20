# MPI 与 Rank 内线程协作分析

本文档重新分析当前项目在 MFEM 并行框架下引入 rank 内线程并行的可行路线。结论先放在前面：

当前项目已经是 MPI 分布式显式矩阵主干，并已加入第一版 rank 内 OpenMP 局部装配执行器。当前路线保持 `ParMesh`、`ParFiniteElementSpace`、`ParBilinearForm` / `ParLinearForm`、`HypreParMatrix` 与现有求解器缓存路径不变，只把本 rank 的局部域单元和边界单元积分、系数求值、局部矩阵/向量计算放到 worker 线程中。

不建议当前阶段直接切到 Partial Assembly、GPU/CEED、MPI 子通信器、场间并发或瞬态力学快照并发。那些路线不是不能做，而是会同时改变矩阵表示、求解器复用方式、输出逻辑和耦合调度，风险远大于收益。

## 1. 当前基线

### 1.1 项目已经完成的 MPI 主干

当前代码中，所有主要求解路径都使用 MFEM 并行对象：

- `main.cpp` 调用 `mfem::Mpi::Init(argc, argv, MPI_THREAD_FUNNELED, &provided)` 与 `mfem::Hypre::Init()`。
- `ConfigLoader::BuildFE()` 读取串行网格后创建 `mfem::ParMesh(MPI_COMM_WORLD, mesh)`。
- 标量场与矢量场分别创建 `mfem::ParFiniteElementSpace`。
- 三个 assembler 使用 `mfem::ParBilinearForm` / `mfem::ParLinearForm`。
- 热场缓存 `mfem::HypreParMatrix` 形式的 `K` 与 `C`，并用 `mfem::Add()` 构造 `K + alpha C`。
- 求解器工厂当前支持 `pcg`、`mumps`、`pardiso`，都走 `HypreParMatrix` 路径。
- 输出阶段通过 `ParMesh::GetSerialMesh()` 和 `ParGridFunction::GetSerialGridFunction()` 在 rank 0 聚合文本结果。

因此，混合并行的第一目标不是“把项目并行化”，而是在现有 MPI-only 并行基础上减少每个 rank 内局部装配和局部系数求值的耗时。

### 1.2 当前已有显式 rank 内装配线程化

当前代码新增了项目自有的 `LocalAssemblyExecutor`，并在静电、热、力学 assembler 中接入。它的并行粒度是本 rank 的域单元和边界单元局部贡献：

- worker 线程各自持有一组 integrator 实例。
- worker 线程并行计算域单元和边界单元上的高斯积分、系数求值和 element matrix/vector。
- 主线程按单元编号顺序做 `DofTransformation` 并把局部贡献合并到 `ParBilinearForm` / `ParLinearForm`。
- 主线程仍独占 MFEM form 的写入、`ParallelAssemble()`、边界消元、线性求解和输出。

因此，本阶段并行的不是全局线性代数求解，也不是 MPI 通信；并行的是当前 MPI rank 内部的域单元与边界单元局部装配计算。

### 1.3 MFEM 提供的相关能力

本仓库中的 MFEM 子模块为 v4.9.1 系列。与本问题直接相关的接口和实现信号有四类。

第一，`mfem::Mpi::Init()` 支持 MPI 线程级别参数。`mfem/general/communication.hpp` 中的接口是：

```cpp
mfem::Mpi::Init(argc, argv, required, &provided);
```

当 `required != MPI_THREAD_SINGLE` 时，MFEM 内部调用 `MPI_Init_thread()`。因此如果项目正式引入 OpenMP worker 线程，入口层不需要绕开 MFEM 手写 MPI 初始化。

第二，MFEM 示例 `examples/ex1p.cpp` 展示了 `Device`、`AssemblyLevel::PARTIAL` 和 `AssemblyLevel::FULL` 的使用方式。这说明 MFEM 原生支持设备后端和部分装配，但这不是当前项目最小改造路线。

第三，`mfem/fem/bilinearform.cpp` 中的 legacy OpenMP 路径只出现在特定元素矩阵缓存函数附近，并且 MFEM 源码明确提示某些 integrator 可能不是线程安全的。这和本项目现状高度相关：自定义表达式、系数和积分子必须先消除共享可变运行时，才能安全进入 rank 内线程装配。

第四，MFEM 的 `FiniteElement` 与 `IntegrationRules` 也需要线程安全构建。当前 `CMakeLists.txt` 在找到 OpenMP 时显式设置 `MFEM_USE_OPENMP=ON`，并始终设置 `MFEM_THREAD_SAFE=ON`。`LocalAssemblyExecutor` 只有在 `_OPENMP`、`MFEM_THREAD_SAFE` 和 `MFEM_USE_OPENMP` 同时可用，且 `OMP_NUM_THREADS > 1` 时才会启用。

### 1.4 当前基线记录

线程安全改造前已记录一个确定性基线：

```bash
env OMP_NUM_THREADS=1 OMP_DYNAMIC=false \
mpirun -np 1 ./build/linux-bundled-release/bin/FEM_solution \
  configs/busbar_etm_nonlinear_transient_vbc.json
```

该基线覆盖 `T` 相关材料表达式、`t` 相关边界表达式、瞬态自适应步长和力学快照批量求解。记录结果：

- accepted steps: `19`
- rejected steps: `5`
- output snapshots: `11`
- baseline file: `results/threadsafe_baseline/busbar_transient_vbc_omp1_np1_prechange.txt`
- `transient.txt` SHA-256: `fd4f3540ad001b2d6b07356dffa3d939dacd37f17cab28c9731ef70fee412654`

当前实现后的验证结果：

- `OMP_NUM_THREADS=1`：走同一本地装配执行器的单 worker 路径，`transient.txt` 与基线 SHA-256 完全一致。
- `OMP_NUM_THREADS=4`：启用 OpenMP 域单元局部装配，accepted/rejected/snapshot 数分别为 `19/5/11`，WRMS 日志与基线一致。
- `OMP_NUM_THREADS=4` 与基线的结果宽表最大绝对差为 `2.0631869710542637e-19`，平均绝对差为 `2.4324804741883336e-21`。
- 边界单元局部计算并行化后，`configs/busbar_electro_thermal_oneway.json` 已在 `OMP_NUM_THREADS=4` 下通过稳态电热运行与 COMSOL 对比。

## 2. 推荐策略

推荐的混合并行策略是：

1. 保留 MPI 分布式显式矩阵主干。
2. 入口层改为请求 `MPI_THREAD_FUNNELED`。
3. 只允许主线程调用 MPI、Hypre、`ParallelAssemble()`、求解器和输出聚合。
4. worker 线程只做本 rank 本地域单元/边界单元的数值积分、局部系数求值和局部贡献生成。
5. 先消除表达式、系数、积分子的共享可变状态，再开启 OpenMP。
6. 线程私有贡献由主线程合并，随后进入 MFEM 的并行装配和线性代数路径。

一句话：MPI 负责把全局问题切开，OpenMP 只负责把每个 rank 自己手里的积分工作算快。

## 3. 不建议优先做的路线

### 3.1 不优先切 Partial Assembly

当前项目依赖显式矩阵的地方很多：

- 热场用 `K`、`C` 两个 `HypreParMatrix` 缓存。
- 自适应步长不断构造 `K + alpha C`。
- Dirichlet 消元使用显式矩阵接口。
- MUMPS/CPARDISO 路径需要显式矩阵。
- 力学批量求解依赖刚度矩阵和求解器缓存。

Partial Assembly 会把核心操作从“显式矩阵 + 直接/迭代求解器”转成“operator action + PA-compatible solver/preconditioner”。这是一条更大的架构迁移，不适合作为当前混合并行第一步。

### 3.2 不优先做场间并发

稳态和瞬态中的 E、T、M 并不是任意独立任务：

- E 可能依赖当前 T 来计算 `sigma(T)`。
- T 依赖 V 和 `sigma` 来计算 Joule 热。
- M 依赖 T 计算热应变。
- Picard 与自适应拒步会改变每个场的有效调用次数。

强行让场间并发会引入调度、同步、状态回滚和内存峰值问题。当前更好的收益点是局部装配热点，而不是场级任务并发。

### 3.3 不优先做瞬态力学快照并发

Phase 2 中不同输出快照的力学 RHS 确实彼此独立，但当前实现复用一个 `MechanicalFieldSolver`、一个温度 `GridFunction` 和一个缓存求解器。并发化会要求：

- 每个线程有独立温度场或只读快照视图。
- 直接求解器是否支持多线程共享因子化需逐一确认。
- 输出快照位移写入需要同步。

这可以作为后续优化，但不是第一阶段。

### 3.4 不请求 `MPI_THREAD_MULTIPLE`

当前目标中 worker 线程不调用 MPI，主线程集中执行通信。因此 `MPI_THREAD_FUNNELED` 足够。`MPI_THREAD_MULTIPLE` 往往成本更高，也会把问题扩大到 MPI 实现、Hypre 和求解器对象的线程安全。

## 4. 线程安全审计

### 4.1 `frontend::Expression`

当前实现已保留 muparser，但把可变求值状态从 `Expression` 中拆出：

- `Expression` 只保存表达式字符串、常量判定和常量值。
- `ExpressionRuntime` 持有 muparser 实例和 `x/y/z/t/T` 变量槽。
- `ExpressionCoefficient` 按线程缓存每个属性对应的 `ExpressionRuntime`。

因此表达式字符串可以共享，parser 与变量槽不能跨线程共享。低频调用仍可直接使用 `Expression::Evaluate()`；积分点热路径应通过线程私有 runtime 缓存求值。

### 4.2 `coeff::ExpressionCoefficient`

当前实现已将只读材料映射与线程私有求值 runtime 分开：

- 共享：属性名、attribute 到表达式字符串的映射、参考温度、温度场指针槽。
- 每线程私有：parser、变量槽、物理坐标缓冲和临时上下文。

温度场通过晚绑定指针读取。进入并行装配区域时，外场必须保持只读，不能由主线程同时更新。

### 4.3 `coeff::JouleHeatingCoefficient`

当前实现已把梯度缓冲改为线程私有。`voltage` 和 `sigma` 可作为只读外场读取，但需要保证线程区域内没有主线程同时修改这些场。

### 4.4 `PiecewiseConstantCoefficient`

该类主要持有 `std::vector<double> attr_to_value_`，构造后只读。它是当前最容易线程安全化的系数类型，只要 `Eval()` 不写共享状态即可。

### 4.5 自定义积分子

以下积分子是 rank 内局部装配的第一批覆盖对象：

- `CustomDiffusionIntegrator`
- `CustomDomainLFIntegrator`
- `CustomBoundaryMassIntegrator`
- `CustomBoundaryLFIntegrator`
- `ThermalStrainLFIntegrator`
- `NormalDisplacementPenaltyIntegrator`

历史实现中的典型共享 scratch 成员包括：

```cpp
mutable mfem::DenseMatrix dshape_buf_, inv_jac_buf_, grad_phys_buf_;
mutable mfem::Vector shape_buf_, normal_buf_;
```

当前已去掉这些成员级 scratch，改为线程私有临时缓冲。OpenMP 装配执行器为每个 worker 使用一份 integrator 实例；这样可以把 MFEM integrator 对象所有权、系数 runtime 和 scratch 生命周期都限制在线程内。

### 4.6 MFEM 对象与通信点

下列操作应留在主线程：

- `ParBilinearForm::ParallelAssemble()`
- `ParLinearForm::ParallelAssemble()`
- `FormLinearSystem()` 中可能触发的并行/约束处理
- `HypreParMatrix` 的全局操作
- `mfem::Add()` 构造并行矩阵
- `CGSolver`、`MUMPSSolver`、`CPardisoSolver`
- `MPI_Allreduce`、`MPI_Barrier`
- ParaView 与 txt 聚合输出

worker 线程区域只生成局部贡献，不进入这些通信或求解阶段。

## 5. 推荐实施路线

### Phase 0：固定运行语义

入口已固定为：

```cpp
int provided = MPI_THREAD_SINGLE;
mfem::Mpi::Init(argc, argv, MPI_THREAD_FUNNELED, &provided);
if (provided < MPI_THREAD_FUNNELED)
{
    spdlog::critical("MPI implementation does not provide MPI_THREAD_FUNNELED");
    return 1;
}
mfem::Hypre::Init();
```

这一步本身不创建线程，只是让 MPI 运行时声明“只有主线程会调用 MPI”。

### Phase 1：增加并行运行配置

建议在 `simulation` 或单独 `parallel` 块中增加：

```jsonc
{
  "parallel": {
    "omp_num_threads": 4,
    "mpi_thread_level": "funneled",
    "rank_local_assembly": false
  }
}
```

当前第一版尚未增加 JSON 并行配置，运行时直接使用 OpenMP 环境变量控制线程数；所有显式装配入口统一走项目自有本地装配执行器。

### Phase 2：表达式运行时线程安全

目标是让材料表达式和时变 BC 表达式都能被多个线程安全求值。

当前做法：

- `Expression` 保留表达式字符串与常量判定。
- 新增 `ExpressionRuntime`，内部持有 parser 和变量槽。
- `ExpressionCoefficient` 每个线程有自己的 runtime cache。
- `UsesVariable("T")` 仍可沿用单线程检测逻辑。

完成后先验证：

- 常量表达式结果不变。
- `electrical_conductivity` 对 `T` 的判断不变。
- 时变电压 BC 在瞬态中仍按 `t` 更新。

### Phase 3：系数线程安全

目标是让 `Coefficient::Eval()` 在 worker 线程里不写共享成员。

当前已完成的基础改造：

1. `ExpressionCoefficient` 的物理坐标缓冲已移到线程私有 runtime cache。
2. `JouleHeatingCoefficient` 的梯度缓冲已改为线程私有。
3. 晚绑定外场指针语义保持不变，并规定并行装配期间外场只读。
4. 只读分域常数系数仍走单独快速路径。

### Phase 4：积分子线程安全

当前已先移除积分子对象中的共享 scratch，并在 OpenMP 装配执行器中采用每线程一份 integrator 实例。这比共享 integrator 再依赖锁更直接，也便于快速证明正确性。

需要覆盖：

- 标量域扩散矩阵。
- 标量域 RHS。
- Robin 边界质量矩阵。
- Robin 边界 RHS。
- 力学热应变 RHS。
- 法向位移罚矩阵。

### Phase 5：本地装配执行器

不要让多个线程同时调用 `SparseMatrix::AddSubMatrix()` 或 MFEM 全局装配容器的可变接口。当前已新增项目自有执行器：

```text
Assembler
  -> LocalAssemblyExecutor
     -> OpenMP parallel for over local attribute/boundary elements
     -> each thread owns integrators, coefficient runtime, scratch buffers
     -> produce per-element matrices/vectors
  -> main thread merges local contributions
  -> ParBilinearForm / ParLinearForm
  -> ParallelAssemble on main thread
```

当前采用第一种落地方式：线程只并行计算域单元和边界单元的 element matrix/vector，主线程按原顺序加入 MFEM form。这样保留 MFEM 现有矩阵、边界消元、`ParallelAssemble()` 和求解器路径，便于与旧结果做严格回归。

## 6. 构建与运行要求

### 6.1 CMake

当前 `CMakeLists.txt` 已查找 OpenMP，把 `OpenMP::OpenMP_CXX` 链接到 `my_fem`，并显式设置：

```cmake
set(MFEM_THREAD_SAFE ON CACHE BOOL "Enable thread-safe MFEM finite element kernels" FORCE)
set(MFEM_USE_OPENMP ON CACHE BOOL "Enable MFEM OpenMP backend" FORCE)
```

这是必要条件：若 MFEM 未以 `MFEM_THREAD_SAFE` 构建，`FiniteElement::CalcShape()` / `CalcDShape()` 可能写共享 scratch；若未启用 `MFEM_USE_OPENMP`，`IntegrationRules` 等全局缓存缺少 OpenMP 保护。当前 release cache 已确认：

```text
MFEM_THREAD_SAFE = YES
MFEM_USE_OPENMP  = YES
```

### 6.2 运行时

推荐测试命令：

```bash
export OMP_NUM_THREADS=4
export OMP_DYNAMIC=false
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpirun -np 2 ./build/linux-bundled-release/bin/FEM_solution configs/busbar_etm_nonlinear_transient_vbc.json
```

资源分配原则：

- `MPI ranks * OMP_NUM_THREADS` 不应超过物理核心数太多。
- 同一个 NUMA 节点内优先让 rank 与其线程贴近。
## 7. 验证路线

正确性验证应先于性能测试。

1. `mpirun -np 1`, `OMP_NUM_THREADS=1`：单 rank、单 worker 路径与当前结果逐项一致。
2. `mpirun -np N`, `OMP_NUM_THREADS=1`：多 rank、单 worker 路径稳定。
3. `mpirun -np N`, `OMP_NUM_THREADS=2/4/8`：同一配置输出一致或在浮点舍入容许范围内一致。
4. 分别验证稳态 E+T、稳态 T+M、瞬态 E+T+M。
5. 对自适应步长场景检查 accepted/rejected step 数、WRMS 日志、输出快照数量是否稳定。
6. 最后再看性能：局部装配耗时、总求解耗时、内存峰值、直接求解器复用是否仍生效。

## 8. 当前不做的事

当前阶段明确不做：

- MPI 子通信器分组。
- 场间任务并发。
- 瞬态 Phase 2 力学快照并发。
- `MPI_THREAD_MULTIPLE`。
- Partial Assembly 主路径迁移。
- GPU / CEED / OCCA 主路径迁移。
- 在线程中调用 MPI、Hypre 或直接求解器。

## 9. 小结

本项目的混合并行路线应从“线程安全的局部装配”开始，而不是从“换 MFEM 后端”或“重排多物理调度”开始。

最重要的工程边界是：

- MPI 通信和求解仍由主线程串行发起。
- worker 线程只做本地数值积分。
- 表达式、系数、积分子必须先去除共享可变运行时。
- 显式矩阵、边界消元、矩阵缓存、直接求解器复用保持不变。

这条路线能最大限度复用当前已经稳定的 MPI/MFEM 主干，同时把风险集中在可测试、可验证的局部计算层。

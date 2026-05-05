# MFEM 常用接口文档

本文档按 MFEM 的常见使用链路整理接口速查，覆盖项目已用和暂未使用的常见对象。MPI/并行专用接口集中放在 [mfem_mpi_interfaces.md](mfem_mpi_interfaces.md)，本项目当前完整使用索引仍保留在 [mfem_api.md](mfem_api.md)。

## 1. 常见对象关系

典型有限元程序的数据链路：

```text
Mesh
  -> FiniteElementCollection
  -> FiniteElementSpace
  -> GridFunction / LinearForm / BilinearForm
  -> SparseMatrix or Operator
  -> Solver
  -> Output
```

| 类别 | 常用类型 | 作用 |
| ---- | -------- | ---- |
| 网格 | `Mesh` | 串行网格、几何和属性信息。 |
| 元集合 | `H1_FECollection`, `L2_FECollection`, `ND_FECollection`, `RT_FECollection` | 定义有限元族和阶次。 |
| 空间 | `FiniteElementSpace` | 将网格和有限元集合组合成自由度空间。 |
| 场变量 | `GridFunction` | 保存节点/自由度上的解。 |
| 系数 | `Coefficient`, `VectorCoefficient` | 为空间积分和边界积分提供物理系数。 |
| 形式 | `BilinearForm`, `LinearForm` | 表达弱形式左端和右端。 |
| 积分子 | `DiffusionIntegrator`, `DomainLFIntegrator` 等 | 将具体弱形式项加入形式对象。 |
| 矩阵向量 | `SparseMatrix`, `Vector`, `Operator` | 线性代数基础对象。 |
| 求解器 | `CGSolver`, `GMRESSolver`, `MINRESSolver`, `DSmoother`, `GSSmoother` | 线性系统求解与预条件。 |
| 输出 | `GridFunction::Save`, `VisItDataCollection`, `ParaViewDataCollection` | 保存结果。 |

## 2. Mesh

### 2.1 构造与读取

```cpp
mfem::Mesh mesh("mesh.msh", 1, 1);
```

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `Dimension()` | 空间维度。 |
| `GetNV()` | 顶点数。 |
| `GetNE()` | 单元数。 |
| `GetNBE()` | 边界单元数。 |
| `GetElement(i)` | 第 `i` 个单元。 |
| `GetBdrElement(i)` | 第 `i` 个边界单元。 |
| `GetElementTransformation(i)` | 单元参考坐标到物理坐标的变换。 |
| `GetElementVertices(i, verts)` | 单元顶点编号。 |
| `GetElementBaseGeometry(i)` | 单元几何类型。 |
| `GetBdrAttribute(i)` | 边界单元属性编号。 |
| `UniformRefinement()` | 全局均匀细化一次。 |
| `Print(stream)` | 输出网格。 |

### 2.2 属性数组

常见成员：

| 成员 | 说明 |
| ---- | ---- |
| `attributes` | 体单元属性集合，通常用于分域材料。 |
| `bdr_attributes` | 边界属性集合，通常用于边界条件。 |

属性编号通常从 `1` 开始；构造 marker 数组时要转换为 `attr - 1` 的下标。

## 3. 有限元集合

常用有限元集合：

| 类型 | 空间 | 常见用途 |
| ---- | ---- | -------- |
| `H1_FECollection(order, dim)` | 连续标量/矢量 H1 | 温度、电势、位移。 |
| `L2_FECollection(order, dim)` | 间断 L2 | 单元常量、DG、投影量。 |
| `ND_FECollection(order, dim)` | Nedelec H(curl) | 电磁场。 |
| `RT_FECollection(order, dim)` | Raviart-Thomas H(div) | 通量、混合方法。 |
| `DG_FECollection(order, dim)` | DG 空间 | 间断 Galerkin。 |

常见模式：

```cpp
auto fec = std::make_unique<mfem::H1_FECollection>(order, mesh.Dimension());
mfem::FiniteElementSpace fes(&mesh, fec.get());
mfem::FiniteElementSpace vfes(&mesh, fec.get(), mesh.Dimension());
```

第三个参数可指定向量维度。向量维度的排序可由可选参数控制，常见值包括 `Ordering::byNODES` 与 `Ordering::byVDIM`。

## 4. FiniteElementSpace

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `GetVSize()` | 本空间向量大小。串行时通常就是自由度数乘向量维度。 |
| `GetNDofs()` | 标量自由度数。 |
| `GetVDim()` | 向量维度。 |
| `GetMesh()` | 关联网格。 |
| `GetFE(i)` | 第 `i` 个单元上的有限元。 |
| `GetElementDofs(i, dofs)` | 获取单元自由度编号。 |
| `GetEssentialTrueDofs(marker, tdofs)` | 根据边界 marker 获取本质边界自由度。串行和并行空间都常用类似接口。 |
| `GetBoundaryTrueDofs(marker, tdofs)` | 获取边界 true dof。 |

使用要点：

- 标量 H1 空间常用于温度、电势。
- 向量 H1 空间常用于位移。
- 对边界条件，先用边界属性构造 marker，再交给空间求本质自由度。

## 5. GridFunction

`GridFunction` 是有限元空间上的场变量。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `GridFunction(fes)` | 在空间上创建场变量。 |
| `ProjectCoefficient(coeff)` | 将系数投影到空间。 |
| `ProjectBdrCoefficient(coeff, marker)` | 在指定边界投影系数。 |
| `GetValue(elem, ip)` | 在单元积分点求值。 |
| `GetVectorValue(elem, ip, vec)` | 在积分点求矢量值。 |
| `GetGradient(elem, ip, grad)` | 标量场梯度。 |
| `SetFromTrueDofs(x)` | 用 true dof 向量更新场变量。 |
| `GetTrueDofs(x)` | 提取 true dof 向量。 |
| `Save(stream)` | 保存场数据。 |

常见用法：

```cpp
mfem::GridFunction u(&fes);
u = 0.0;
u.ProjectCoefficient(coef);
```

## 6. Coefficient

MFEM 通过 coefficient 对象把物理参数传入积分子。

### 6.1 标量系数

| 类型 | 说明 |
| ---- | ---- |
| `ConstantCoefficient` | 常数标量。 |
| `PWConstCoefficient` | 按属性分片常数。 |
| `FunctionCoefficient` | C++ 函数回调标量。 |
| `GridFunctionCoefficient` | 从已有场变量取值。 |
| `ProductCoefficient` | 两个标量系数乘积。 |
| `SumCoefficient` | 两个标量系数求和。 |
| `TransformedCoefficient` | 对系数做函数变换。 |

示例：

```cpp
mfem::ConstantCoefficient k(400.0);
mfem::FunctionCoefficient src([](const mfem::Vector &x) { return 1.0; });
```

### 6.2 矢量和矩阵系数

| 类型 | 说明 |
| ---- | ---- |
| `VectorConstantCoefficient` | 常矢量系数。 |
| `VectorFunctionCoefficient` | 函数回调矢量。 |
| `VectorGridFunctionCoefficient` | 从矢量场取值。 |
| `MatrixConstantCoefficient` | 常矩阵系数。 |
| `MatrixFunctionCoefficient` | 函数回调矩阵。 |

自定义系数通常继承 `Coefficient` 或 `VectorCoefficient` 并重写 `Eval()`。

## 7. BilinearForm

`BilinearForm` 表示弱形式左端。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `AddDomainIntegrator(integ)` | 添加体积分项。 |
| `AddBoundaryIntegrator(integ)` | 添加边界积分项。 |
| `AddInteriorFaceIntegrator(integ)` | 添加内部面项，常用于 DG。 |
| `Assemble()` | 装配。 |
| `Finalize()` | 完成稀疏矩阵构造。 |
| `SpMat()` | 获取串行稀疏矩阵。 |
| `FormLinearSystem(ess_tdof, x, b, A, X, B)` | 根据边界条件生成线性系统。 |
| `RecoverFEMSolution(X, b, x)` | 从线性系统解恢复有限元解。 |

常见积分子：

| 积分子 | 弱形式项 |
| ------ | -------- |
| `DiffusionIntegrator(k)` | `k grad(u) . grad(v)` |
| `MassIntegrator(c)` | `c u v` |
| `ConvectionIntegrator(a)` | 对流项。 |
| `VectorDiffusionIntegrator(c)` | 矢量扩散/弹性近似项。 |
| `ElasticityIntegrator(lambda, mu)` | 线弹性刚度。 |
| `BoundaryMassIntegrator(h)` | 边界质量项，如 Robin 左端。 |

## 8. LinearForm

`LinearForm` 表示弱形式右端。

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `AddDomainIntegrator(integ)` | 添加体源项。 |
| `AddBoundaryIntegrator(integ)` | 添加边界源项。 |
| `Assemble()` | 装配右端向量。 |
| `Update()` | 空间变化后更新内部数据。 |

常见积分子：

| 积分子 | 弱形式项 |
| ------ | -------- |
| `DomainLFIntegrator(f)` | `f v` |
| `BoundaryLFIntegrator(g)` | `g v` |
| `VectorDomainLFIntegrator(f)` | 矢量体力。 |
| `VectorBoundaryLFIntegrator(g)` | 矢量边界力。 |

## 9. 矩阵、向量和 Operator

### 9.1 Vector

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `SetSize(n)` | 设置长度。 |
| `Size()` | 当前长度。 |
| `operator=(value)` | 全部赋值。 |
| `Norml2()` | L2 范数。 |
| `Max()` / `Min()` | 最大/最小值。 |
| `Add(a, y)` | 执行 `x += a y`。 |
| `Print(stream)` | 输出。 |

### 9.2 SparseMatrix

| 接口 | 说明 |
| ---- | ---- |
| `Height()` / `Width()` | 矩阵尺寸。 |
| `Mult(x, y)` | 计算 `y = A x`。 |
| `Add(a, B)` | 加入另一矩阵。 |
| `Finalize()` | 完成装配。 |

### 9.3 Operator

`Operator` 是 MFEM 求解器和矩阵的共同抽象。

| 接口 | 说明 |
| ---- | ---- |
| `Height()` / `Width()` | 算子尺寸。 |
| `Mult(x, y)` | 算子乘法。 |
| `SetOperator(op)` | 求解器设置待求解算子。 |

很多矩阵、预条件器和求解器都可作为 `Operator` 使用。

## 10. 求解器和预条件器

常用迭代求解器：

| 类型 | 适用场景 |
| ---- | -------- |
| `CGSolver` | 对称正定系统。 |
| `GMRESSolver` | 一般非对称系统。 |
| `MINRESSolver` | 对称不定系统。 |
| `BiCGSTABSolver` | 一般非对称系统。 |

常用设置：

```cpp
mfem::CGSolver cg;
cg.SetRelTol(1e-12);
cg.SetAbsTol(0.0);
cg.SetMaxIter(500);
cg.SetPrintLevel(1);
cg.SetOperator(A);
cg.Mult(b, x);
```

常用预条件器：

| 类型 | 说明 |
| ---- | ---- |
| `DSmoother` | 对角/Jacobi 平滑。 |
| `GSSmoother` | Gauss-Seidel 平滑。 |
| `OperatorJacobiSmoother` | 基于 operator 的 Jacobi。 |
| `BlockDiagonalPreconditioner` | 分块对角预条件。 |
| `BlockLowerTriangularPreconditioner` | 分块下三角预条件。 |

## 11. 时间相关接口

MFEM 提供 ODE/time-dependent 框架，适合把半离散系统写成 `du/dt = f(u,t)`。

常见类型：

| 类型 | 说明 |
| ---- | ---- |
| `TimeDependentOperator` | 时间相关算子基类。 |
| `ODESolver` | ODE 求解器基类。 |
| `ForwardEulerSolver` | 显式 Euler。 |
| `RK2Solver`, `RK3SSPSolver`, `RK4Solver` | Runge-Kutta 系列。 |
| `BackwardEulerSolver` | 隐式 Euler。 |
| `SDIRK23Solver`, `SDIRK33Solver` | 隐式 Runge-Kutta。 |

本项目目前手写热方程 BDF 时间推进，没有直接使用 MFEM 的 ODE 框架。

## 12. 输出接口

### 12.1 简单保存

```cpp
std::ofstream mesh_ofs("mesh.mesh");
mesh.Print(mesh_ofs);

std::ofstream sol_ofs("sol.gf");
u.Save(sol_ofs);
```

### 12.2 DataCollection

常用类型：

| 类型 | 说明 |
| ---- | ---- |
| `VisItDataCollection` | VisIt/GLVis 风格输出。 |
| `ParaViewDataCollection` | ParaView `.pvd/.vtu` 输出。 |

常用接口：

| 接口 | 说明 |
| ---- | ---- |
| `RegisterField(name, gf)` | 注册场变量。 |
| `SetPrefixPath(path)` | 设置输出目录。 |
| `SetCycle(cycle)` | 设置步编号。 |
| `SetTime(time)` | 设置物理时间。 |
| `Save()` | 写入文件。 |

## 13. 设备与装配等级

MFEM 支持不同计算后端和装配策略。

### 13.1 Device

```cpp
mfem::Device device("cpu");
device.Print();
```

常见后端字符串包括 `cpu`、`cuda`、`hip`、`omp` 等，具体可用性取决于 MFEM 构建选项。

### 13.2 AssemblyLevel

| 枚举 | 说明 |
| ---- | ---- |
| `AssemblyLevel::LEGACY` | 传统全矩阵装配。 |
| `AssemblyLevel::FULL` | 完整装配。 |
| `AssemblyLevel::PARTIAL` | Partial Assembly，常用于 GPU/高阶。 |
| `AssemblyLevel::ELEMENT` | 元级装配。 |
| `AssemblyLevel::NONE` | 矩阵自由路径。 |

本项目当前主路径是显式矩阵装配，没有把 Partial Assembly 作为默认架构。

## 14. 调试与检查

常用辅助接口：

| 接口 | 说明 |
| ---- | ---- |
| `MFEM_VERIFY(cond, msg)` | 条件检查，失败时报错。 |
| `MFEM_ABORT(msg)` | 直接中止并输出错误。 |
| `out` / `err` | MFEM 输出流封装。 |
| `Memory<T>` / `Array<T>` | MFEM 内存与数组工具。 |

## 15. 选择建议

| 需求 | 建议接口 |
| ---- | -------- |
| 连续标量扩散问题 | `H1_FECollection` + `BilinearForm` + `DiffusionIntegrator`。 |
| 热瞬态质量矩阵 | `MassIntegrator`。 |
| 线弹性 | `ElasticityIntegrator` 或自定义弹性积分子。 |
| 有分域材料 | `PWConstCoefficient` 或自定义按属性系数。 |
| 对已有场求源项 | `GridFunctionCoefficient` 或自定义系数。 |
| 串行快速原型 | `Mesh` + `FiniteElementSpace` + `SparseMatrix`。 |
| MPI 生产路径 | 见 [mfem_mpi_interfaces.md](mfem_mpi_interfaces.md)。 |

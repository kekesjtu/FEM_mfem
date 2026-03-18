<!-- markdownlint-disable MD024 -->
# FEM_mfem 框架接口文档（仅覆盖自研层）

## 1. 文档范围

本文档描述当前仓库自研框架层（`include/fem/*` 与 `src/*`）的接口与行为约定，不展开 MFEM 子仓库内部 API。

- 覆盖：配置解析、表达式系统、材料库、FE 上下文、组装、求解、后处理、应用编排
- 不覆盖：`mfem/` 目录内部实现

> 文档基线：按当前主工程代码状态更新（截至 2026-03-18）。

---

## 2. 分层与依赖规则

### 2.1 依赖方向

- `frontend` / `material`：只表达“物理定义与参数”，不依赖具体求解实现。
- `fe` / `mapping` / `assembly` / `solver`：后端执行层，不依赖 JSON 解析细节。
- `app`：编排层，负责按场类型选择流水线并串联前后端。
- `post` / `io`：求解后结果导出与后处理，不反向依赖 `frontend`。

### 2.2 命名空间

- `fem::frontend`：配置模型、配置加载、表达式求值
- `fem::material`：材料属性数据库
- `fem::fe`：有限元上下文（标量/力学）
- `fem::mapping`：表达式到 MFEM 系数映射
- `fem::assembly`：组装抽象与 MFEM 适配器
- `fem::solver`：线性求解器抽象与实现
- `fem::post`：文本导出与应力后处理
- `fem::io`：ParaView 导出
- `fem::log`：日志初始化与访问
- `fem::app`：应用入口、场调度、耦合流水线

---

## 3. 前端层接口（frontend）

## 3.1 `SimulationConfig`

职责：全局求解控制。

- `mesh_path`：网格路径（必填）
- `order = 1`：有限元阶次
- `uniform_refinement_levels = 0`：统一加密次数
- `output_dir = "results"`：输出目录
- `log_level = "info"`：日志级别

## 3.2 `MaterialConfig`

职责：材料定义。

- `name`：材料名
- `properties`：属性名 -> 表达式字符串

常见属性键（按流程使用）：

- 标量扩散类：`diffusion`
- 电-热：`electrical_conductivity`（缺失时可回退 `diffusion`）
- 力学：`young_modulus`、`poisson_ratio`
- 热应变：`thermal_expansion_secant`（兼容旧键 `thermal_expansion`）

## 3.3 `FieldConfig`

职责：单个物理场定义；支持 `electrostatic` / `thermal` / `mechanical`。

标量场相关字段：

- `source` / `source_default` / `source_by_domain`
- `dirichlet_value` / `dirichlet_bdr_attributes`
- `dirichlet_boundary_conditions`（分边界常值 Dirichlet）
- `robin_l` / `robin_q` / `robin_bdr_attributes`
- `robin_on_remaining_boundaries`
- `robin_boundary_conditions`（分边界 Robin）

力学场相关字段：

- `body_force_default` / `body_force_by_domain`
- `dirichlet_displacement_value`
- `traction_boundary_conditions`
- `pressure_boundary_conditions`
- `normal_displacement_boundary_conditions`（罚项弱施加 $u\cdot n = 0$）
- `reference_temperature`
- `enable_thermal_strain_coupling`
- `enable_stress_postprocess`

## 3.4 `ProjectConfig`

- `simulation`
- `materials`
- `domain_materials`（材料名 -> 域属性号列表）
- `fields`

## 3.5 `ConfigLoader`

接口：

- `static ProjectConfig LoadFromFile(const std::string &path);`

行为：

- 从 JSON 读入 `ProjectConfig`。
- 按默认值补齐可选字段。
- 读取 `domain_materials`、`source_by_domain`、各种边界条件与力学扩展字段。

异常：

- 文件不可读：`std::runtime_error("Failed to open config file: ...")`
- 必填键缺失（如 `simulation.mesh_path`、`fields[*].name/type`）：JSON 访问异常上抛。

## 3.6 `EvalContext` 与 `ScalarExpression`

- `EvalContext`：`x,y,z,t` + `coupled_fields`
- `ScalarExpression`：支持常量快速路径与 muparser 通用表达式。
- 构造时可注册耦合变量名；求值时未提供值默认按 `0.0`。
- `Raw()` 返回原始表达式文本。

线程安全：单个 `ScalarExpression` 内含可变缓存，不保证并发线程安全。

---

## 4. 材料层接口（material）

## 4.1 `MaterialDatabase`

接口：

- `Add(const frontend::MaterialConfig&, const std::vector<std::string>& coupled_variable_names = {})`
- `GetProperty(const std::string &material, const std::string &property) const`
- `Contains(const std::string &material, const std::string &property) const`

行为：

- `Add` 时将属性表达式预编译为 `ScalarExpression`。
- `GetProperty` 返回只读引用；不存在则抛异常。

异常：

- 材料不存在：`Material not found: ...`
- 属性不存在：`Property not found: material.property`

---

## 5. FE 基础层接口（fe）

## 5.1 `ScalarFeContext`

接口：

- `ScalarFeContext(const std::string &mesh_path, int order, int uniform_refinement_levels)`
- `mfem::Mesh &Mesh()`
- `mfem::FiniteElementSpace &Space()`

行为：

- 加载网格、执行统一加密、构建标量 `H1` 空间。

## 5.2 `MechanicalFeContext`

接口：

- `MechanicalFeContext(...)`
- `mfem::Mesh &Mesh()`
- `mfem::FiniteElementSpace &Space()`
- `int Dimension() const`

行为：

- 构建向量位移场对应的 `H1` 空间（`vdim = mesh.Dimension()`）。

生命周期约束：`Mesh` 与 `FESpace` 引用不可越过上下文对象生命周期。

---

## 6. 映射层接口（mapping）

## 6.1 `CoefficientFactory`

接口：

- `CreateScalar(const fem::frontend::ScalarExpression &expr)`

行为：

- 生成 `mfem::Coefficient`（`FunctionCoefficient` 适配），将积分点坐标注入 `EvalContext.x/y/z`。

---

## 7. 组装层接口（assembly）

## 7.1 Poisson 组装抽象

核心类型：

- `DirichletBoundaryCondition`（常值 + marker）
- `RobinBoundaryCondition`（`l,q` + marker）
- `PoissonAssemblyInput`
- `AssembledSystem`
- `IPoissonAssembler`
- `MfemPoissonAssembler`

`MfemPoissonAssembler` 行为：

- 双线性项：`CustomDiffusionIntegrator` + 可选 Robin 质量项
- 线性项：`CustomDomainLFIntegrator` + 可选 Robin 载荷项
- 计算本质边界自由度后 `FormLinearSystem`
- 若提供分边界 Dirichlet 条件，按 marker 投影初值；否则使用 `dirichlet_value`

## 7.2 线弹性组装抽象

核心类型：

- `MechanicalTractionBoundary`
- `MechanicalPressureBoundary`
- `MechanicalNormalDisplacementBoundary`
- `LinearElasticityInput`
- `MechanicalSystem`
- `IMechanicalAssembler`
- `MfemLinearElasticityAssembler`

`MfemLinearElasticityAssembler` 行为：

- 双线性项：`ElasticityIntegrator` + 可选法向位移罚项
- 线性项：体力 + 牵引 + 压力 + 可选热应变载荷 `ThermalStrainLFIntegrator`
- 按位移 Dirichlet 投影边界初值并形成线性系统

## 7.3 自定义积分器

- `PoissonIntegrators.hpp`：标量扩散/源项/边界项自定义积分器
- `LinearElasticityIntegrators.hpp`：
 	- `ThermalStrainLFIntegrator`
 	- `NormalDisplacementPenaltyIntegrator`

---

## 8. 求解层接口（solver）

## 8.1 `ILinearSolver`

- `virtual void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) = 0;`

## 8.2 `MfemPcgSolver`

构造：

- `MfemPcgSolver(double rel_tol, double abs_tol, int max_iter, int print_level)`

行为：

- 要求 `A` 可转换为 `mfem::SparseMatrix`
- 预条件：`mfem::DSmoother`
- 迭代：`mfem::PCG`

异常：

- `A` 非 `SparseMatrix` 时抛 `std::runtime_error`

---

## 9. 输出与后处理接口

## 9.1 `fem::io::ParaviewExporter`

接口：

- `ParaviewExporter(std::string output_dir)`
- `Save(collection_name, field_name, mesh, field, cycle, time)`

行为：

- 自动创建输出目录，使用 `mfem::ParaViewDataCollection` 输出。

## 9.2 `fem::post::SolutionTextExporter`

接口：

- `ExportScalarNodalTxt(...)`
- `ExportVectorNodalTxt(...)`

行为：

- 按网格顶点采样写出 txt。
- 向量导出包含分量与模长。

## 9.3 `fem::post::MechanicalPostProcessor`

接口：

- `FillVonMisesNodalField(...)`
- `ExportElementStressCsv(...)`

行为：

- 计算节点 von Mises（数值积分加权）。
- 导出单元中心应力 CSV（`sxx/syy/szz/sxy/syz/szx/von_mises`）。
- 可选热应变修正（通过 `ThermalStrainParameters`）。

---

## 10. 日志接口（log）

`namespace fem::log`：

- `Init(const std::string &level)`
- `Get()`

支持级别：`trace/debug/info/warn/error/critical`。

---

## 11. 应用编排层接口（app）

## 11.1 入口

- `class Application { int Run(const std::string &config_path); }`
- `main.cpp` 默认配置路径：`configs/bipolar_thermo_mechanical.json`

## 11.2 Pipeline 选择与调度

`ApplicationPipelines.hpp` 核心能力：

- 提取耦合变量名：`BuildCoupledVariableNames`
- 材料库构建：`BuildMaterialDatabase`
- 物理场筛选：`SelectFieldsByType`
- 优先尝试耦合流水线：`RunCoupledPipelineByPhysicsType`
- 否则按场类型逐个执行：`RunSingleFieldPipelineByPhysicsType`

耦合策略（当前实现）：

- 若存在 `thermal + mechanical`：执行热-力单向耦合
- 若存在 `electrostatic + thermal` 且无力学场：执行电-热单向耦合
- 其余情况：按单场 runner 逐场执行

## 11.3 Field Runner 体系

- `CreateFieldRunners()` 映射：
 	- `electrostatic` -> `RunScalarField`
 	- `thermal` -> `RunScalarField`
 	- `mechanical` -> `RunMechanicalField`

### 标量场流程（`RunScalarField`）

1. 构建 `ScalarFeContext`
2. 校验 `domain_materials` 覆盖所有网格域
3. 构建分段扩散系数（来自材料 `diffusion`）
4. 构建分段源项系数（`source_default/source_by_domain`）
5. 建边界条件并调用 Poisson 组装 + PCG 求解
6. 输出 txt 与 ParaView，并记录统计量

### 力学场流程（`RunMechanicalFieldInternal`）

1. 构建 `MechanicalFeContext`
2. 校验域覆盖
3. 从材料 `young_modulus/poisson_ratio` 构建分段 Lamé 系数
4. 组装体力、位移边界、牵引/压力/法向位移罚项
5. 可选热应变耦合（`alpha + temperature`）
6. PCG 求解位移场
7. 可选 von Mises 与单元应力 CSV 后处理导出

### 电-热单向耦合（`RunElectroThermalCoupled`）

1. 先解电势场
2. 计算焦耳热项 $Q=\sigma\|\nabla V\|^2$
3. 与热源叠加后再解温度场
4. 导出双场结果

### 热-力单向耦合（`RunThermoMechanicalCoupled`）

1. 先解温度场
2. 温度场投影为力学热应变输入
3. 解位移场并执行力学后处理

---

## 12. 配置契约（实现对应）

对应 schema：`docs/physics.schema.json`。

关键约束（按当前代码）：

- `simulation.mesh_path` 必填
- `fields` 至少 1 项，且 `fields[*].type ∈ {electrostatic, thermal, mechanical}`
- 实际运行依赖 `domain_materials` 覆盖全部网格域（标量/耦合/力学路径均会校验）
- 标量扩散系数来自材料属性 `diffusion`
- 电导率优先 `electrical_conductivity`，缺失时回退 `diffusion`
- 力学材料需提供 `young_modulus` 与 `poisson_ratio`
- 热应变系数优先 `thermal_expansion_secant`，兼容旧键 `thermal_expansion`

---

## 13. 当前能力边界

- 当前耦合为单向顺序耦合，未实现迭代收敛闭环。
- 暂无瞬态时间推进器（仅稳态流程）。
- `ScalarExpression` 实例级别非线程安全。
- `domain_materials` 在运行期是强约束：未覆盖域会抛异常终止。

---

## 14. 版本说明

当前文档对应代码能力：

- 标量 Poisson（电/热）统一流程
- 3D 线弹性（体力、牵引、压力、法向位移罚项）
- 电-热单向耦合与热-力单向耦合
- muparser 表达式系统（含耦合变量）
- 分域材料映射、Robin 边界、解统计日志
- 力学 von Mises 后处理与单元应力 CSV 导出

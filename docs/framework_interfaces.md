<!-- markdownlint-disable MD024 -->
# FEM_mfem 框架接口文档（仅覆盖自研层）

## 1. 文档范围

本文件仅说明当前工程中新建的框架层接口与实现约定，不展开 MFEM 原生 API 细节。

- 覆盖目录：`include/fem/*` 与 `src/*`
- 不覆盖：`mfem/` 子仓库内部接口

---

## 2. 分层与依赖规则

### 2.1 依赖方向

- 前端层（`frontend`, `material`）只表达“物理定义与参数”，不依赖求解实现。
- 后端层（`assembly`, `solver`）只依赖抽象输入与 MFEM 适配器，不依赖 JSON 解析细节。
- `app` 作为编排层，负责把前端定义映射到后端对象。

### 2.2 命名空间

- `fem::frontend`：配置与表达式
- `fem::material`：材料库
- `fem::fe`：有限元基础上下文（Mesh/FESpace）
- `fem::mapping`：前端表达式到 MFEM Coefficient 映射
- `fem::assembly`：组装抽象与 MFEM 适配器
- `fem::solver`：求解器抽象与 MFEM 适配器
- `fem::io`：输出
- `fem::log`：日志
- `fem::app`：应用编排入口

---

## 3. 前端层接口（frontend）

## 3.1 `struct fem::frontend::SimulationConfig`

**职责**：保存全局求解控制参数。

字段：

- `std::string mesh_path`：网格文件路径（必填）
- `int order = 1`：有限元阶次
- `int uniform_refinement_levels = 0`：全局加密次数
- `std::string output_dir = "results"`：输出目录
- `std::string log_level = "info"`：日志级别

## 3.2 `struct fem::frontend::MaterialConfig`

**职责**：描述一个材料及其属性表达式。

字段：

- `std::string name`：材料名
- `std::unordered_map<std::string, std::string> properties`：属性名 -> 表达式字符串

约定：

- 当前 Poisson 型流程默认读取属性键 `diffusion`。

## 3.3 `struct fem::frontend::FieldConfig`

**职责**：描述一个物理场实例。

字段：

- `std::string name`：场变量名（可被其他表达式耦合引用）
- `std::string type`：场类型（当前支持 `electrostatic` / `thermal` / `mechanical`）
- `std::string material`：绑定材料名
- `std::string source = "0.0"`：源项表达式
- `std::string source_default = "0.0"`：默认体源项（域分段时的默认值）
- `std::unordered_map<int, std::string> source_by_domain`：按域属性号覆盖源项
- `std::string diffusion = "1.0"`：默认扩散/导热系数表达式
- `std::unordered_map<int, std::string> diffusion_by_domain`：按域属性号覆盖扩散系数
- `double dirichlet_value = 0.0`：Dirichlet 常值
- `std::vector<int> dirichlet_bdr_attributes`：边界属性号（1-based）
- `double robin_l = 0.0`：Robin 边界左端系数 $L$
- `double robin_q = 0.0`：Robin 边界右端项 $q$
- `std::vector<int> robin_bdr_attributes`：Robin 显式边界属性号（1-based）
- `bool robin_on_remaining_boundaries = false`：Robin 是否施加到除 Dirichlet 外的其余边界

## 3.4 `struct fem::frontend::ProjectConfig`

**职责**：项目总配置根对象。

字段：

- `SimulationConfig simulation`
- `std::vector<MaterialConfig> materials`
- `std::unordered_map<std::string, std::vector<int>> domain_materials`：材料到域属性号列表映射
- `std::vector<FieldConfig> fields`

## 3.5 `class fem::frontend::ConfigLoader`

### 接口

- `static ProjectConfig LoadFromFile(const std::string &path);`

### 行为

- 从 JSON 文件读取配置并填充 `ProjectConfig`。
- 对缺省项写入默认值（如 `order`, `output_dir`, `source` 等）。
- 支持读取 `domain_materials`、`source_by_domain`、`diffusion_by_domain` 与 Robin 参数。

### 异常

- 文件不可读：抛出 `std::runtime_error`
- 缺少关键字段（例如 `simulation.mesh_path`）：JSON 访问时抛异常

## 3.6 `struct fem::frontend::EvalContext`

**职责**：表达式求值上下文。

字段：

- `double x,y,z,t`：时空变量
- `std::unordered_map<std::string,double> coupled_fields`：耦合场变量值

## 3.7 `class fem::frontend::ScalarExpression`

### 设计目标

- 统一表达式载体，支持常数优化 + muparser 通用表达式。
- 支持内置变量 `x,y,z,t` 与耦合变量（字段名）。

### 构造与生命周期

- `explicit ScalarExpression(std::string expression, std::vector<std::string> coupled_variable_names = {});`
- 支持复制、移动、析构（已显式定义）。

### 核心方法

- `double Evaluate(const EvalContext &ctx) const;`
- `const std::string &Raw() const;`

### 行为约定

- 若表达式是纯数值常量，走常量快速路径。
- 否则通过 muparser 解析；未提供的耦合变量在求值时默认视为 `0.0`。
- muparser 错误统一转为 `std::runtime_error`。

### 线程安全说明

- 当前对象内部含可变状态（变量缓存），**同一个 `ScalarExpression` 实例不保证并发线程安全**。
- 并行场景建议：线程内独立实例，或上层加锁。

---

## 4. 材料层接口（material）

## 4.1 `class fem::material::MaterialDatabase`

### 接口

- `void Add(const frontend::MaterialConfig &config, const std::vector<std::string> &coupled_variable_names = {});`
- `const frontend::ScalarExpression &GetProperty(const std::string &material, const std::string &property) const;`

### 行为

- `Add` 将材料属性字符串预编译为 `ScalarExpression`。
- `GetProperty` 返回属性表达式只读引用。

### 异常

- 材料不存在：`std::runtime_error("Material not found: ...")`
- 属性不存在：`std::runtime_error("Property not found: ...")`

---

## 5. FE 基础层接口（fe）

## 5.1 `class fem::fe::FeContext`

### 接口

- `FeContext(const std::string &mesh_path, int order, int uniform_refinement_levels);`
- `mfem::Mesh &Mesh();`
- `mfem::FiniteElementSpace &Space();`

### 行为

- 负责网格加载、统一加密、`H1_FECollection` 与 `FiniteElementSpace` 构造。
- 遵循“最小封装”原则，直接暴露 MFEM 原生对象引用。

### 生命周期约束

- `Mesh` 与 `FiniteElementSpace` 的生命周期由 `FeContext` 管理。
- 任何外部引用不得超过 `FeContext` 生命周期。

---

## 6. 映射层接口（mapping）

## 6.1 `class fem::mapping::CoefficientFactory`

### 接口

- `static std::unique_ptr<mfem::Coefficient> CreateScalar(const fem::frontend::ScalarExpression &expr);`

### 行为

- 将 `ScalarExpression` 映射为 MFEM `FunctionCoefficient`。
- 自动把空间坐标映射到 `EvalContext.x/y/z`。

### 扩展点

- 可新增：向量系数、矩阵系数工厂方法。
- 可新增：依赖 `t` 与多场状态的动态回调输入。

---

## 7. 组装层接口（assembly）

## 7.1 `struct fem::assembly::PoissonAssemblyInput`

字段：

- `mfem::FiniteElementSpace &space`
- `mfem::Coefficient &diffusion`
- `mfem::Coefficient &source`
- `mfem::Array<int> essential_bdr_marker`
- `mfem::Coefficient *robin_l = nullptr`
- `mfem::Coefficient *robin_q = nullptr`
- `mfem::Array<int> robin_bdr_marker`
- `double dirichlet_value = 0.0`

## 7.2 `struct fem::assembly::AssembledSystem`

字段：

- `std::unique_ptr<mfem::BilinearForm> bilinear`
- `std::unique_ptr<mfem::LinearForm> linear`
- `mfem::OperatorPtr A`
- `mfem::Vector X`
- `mfem::Vector B`
- `mfem::Array<int> essential_tdofs`

## 7.3 `class fem::assembly::IAssembler`

### 接口

- `virtual AssembledSystem Assemble(PoissonAssemblyInput &, mfem::GridFunction &initial_guess) = 0;`

### 设计意图

- 抽象组装行为，隔离“方程定义”与“具体 MFEM 组装实现”。
- 后续可替换为自研内核，不改上层业务流程。

## 7.4 `class fem::assembly::MfemPoissonAssembler final : public IAssembler`

### 行为

- 组装 `DiffusionIntegrator` 与 `DomainLFIntegrator`。
- 当提供 Robin 参数时，额外组装边界质量项与边界载荷项（`MassIntegrator` / `BoundaryLFIntegrator`）。
- 计算本质边界自由度，调用 `FormLinearSystem` 形成 `A X = B`。

---

## 8. 求解层接口（solver）

## 8.1 `class fem::solver::ILinearSolver`

### 接口

- `virtual void Solve(mfem::Operator &A, const mfem::Vector &b, mfem::Vector &x) = 0;`

### 设计意图

- 解耦线性系统解法与上层流程。
- 允许在不动组装层的前提下切换求解器。

## 8.2 `class fem::solver::MfemPcgSolver final : public ILinearSolver`

### 构造

- `MfemPcgSolver(double rel_tol, double abs_tol, int max_iter, int print_level);`

### 行为

- 当前实现假设 `A` 可转为 `mfem::SparseMatrix`。
- 默认使用 `mfem::GSSmoother` 预条件 + `mfem::PCG`。

### 异常

- `A` 不是 `SparseMatrix`：抛 `std::runtime_error`

---

## 9. 输出层接口（io）

## 9.1 `class fem::io::ParaviewExporter`

### 接口

- `explicit ParaviewExporter(std::string output_dir);`
- `void Save(const std::string &collection_name, const std::string &field_name, mfem::Mesh &mesh, mfem::GridFunction &field, int cycle, double time) const;`

### 行为

- 直接使用 `mfem::ParaViewDataCollection` 输出。
- 自动创建输出目录。

---

## 10. 日志接口（log）

## 10.1 `namespace fem::log`

### 接口

- `void Init(const std::string &level);`
- `std::shared_ptr<spdlog::logger> Get();`

### 行为

- 初始化默认彩色控制台 logger。
- 支持级别：`trace/debug/info/warn/error/critical`。

---

## 11. 应用编排接口（app）

## 11.1 `class fem::app::Application`

### 接口

- `int Run(const std::string &config_path);`

### 当前流程（首版）

1. 加载 JSON 配置
2. 初始化日志
3. 收集场名并注册为耦合变量名
4. 构建材料库与 FE 上下文
5. 取首个场并检查类型
6. 构造默认/分段扩散系数与源项系数
7. 构造 Dirichlet 与 Robin 边界标记
8. 组装线性系统（含可选 Robin 边界项）
9. 调用线性求解器并恢复网格解
10. 计算场统计量（min/max/mean/non-finite）并导出 ParaView

### 当前限制

- 只执行 `fields.front()`（即单场入口）
- 多场耦合目前仍是“表达式变量层面可引用”，尚未形成迭代/分块求解闭环
- 暂未实现瞬态时间推进器

---

## 12. 配置文件契约（实现对应）

对应 schema：`docs/physics.schema.json`

最小示例：`configs/steady_electrostatic.json`

关键约束：

- `simulation.mesh_path` 必填
- `fields` 至少 1 项
- `fields[*].type` 需在支持列表内
- 若设置 `field.material`，则材料中需存在相应属性（当前要求 `diffusion`）
- 热问题可选使用：`domain_materials` + `source_by_domain` + Robin 参数

---

## 13. 扩展建议（保持接口稳定）

### 13.1 多场耦合求解

- 在 `app` 层增加 `IFieldCouplingStrategy`（分块迭代、牛顿耦合等）
- 将 `fields` 循环与耦合变量回写机制抽象化

### 13.2 瞬态支持

- 新增 `ITimeIntegrator` 抽象（BE/CN/RK）
- `assembly` 层补质量矩阵接口，支持 `M dU/dt + K U = F`

### 13.3 非线性支持

- 新增 `INonlinearAssembler` 与 `INonlinearSolver` 接口
- 复用现有前端表达式系统作为系数回调来源

---

## 14. 版本说明

- 文档对应代码状态：DIP 架构骨架 + Poisson/稳态热传导统一工作流 + muparser 表达式 + 分段域系数 + Robin 边界 + 解统计日志

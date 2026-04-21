# 自适应时间步长策略详细说明

## 1. 总体架构

本求解器采用 **预测-校正（Predictor-Corrector）** 框架，结合 **自适应阶数选择** 和 **WRMS（Weighted Root Mean Square）误差控制**，实现瞬态热传导问题的自适应时间积分。

核心思想：每一步同时计算两个解——预测值（外推）和校正值（隐式求解），二者之差作为局部截断误差的估计，据此动态调整下一步的时间步长 $\Delta t$。

### 时间积分方法

| 阶数 | 预测器（Predictor） | 校正器（Corrector） | 适用阶段            |
| ---- | ------------------- | ------------------- | ------------------- |
| 1 阶 | 线性外推            | 后向欧拉（BE）      | `history_depth ≥ 1` |
| 2 阶 | 二次 Lagrange 外推  | 变步长 BDF2         | `history_depth ≥ 2` |

---

## 2. 冷启动与一致性初始化

### 2.1 初始时间步长选取

采用 Hairer-Wanner 风格的初始步长估计：

1. 计算初始变化率：

$$\dot{T}_0 \approx \text{diag}(C)^{-1}(F_0 - K \cdot T_0)$$

取全局最大值 `rate = max |dT/dt|₀`。

2. 估计特征尺度：

$$d_0 = \text{abstol} + \text{reltol} \cdot |T_{\text{init}}|$$

$$d_1 = \max\left(\text{rate},\; \frac{d_0}{\Delta t_{\max}}\right)$$

3. 初始步长：

$$\Delta t_{\text{init}} = \text{clamp}\left(0.01 \cdot \frac{d_0}{d_1},\; \Delta t_{\min},\; \Delta t_{\max}\right)$$

**物理含义**：
- 当 `rate ≈ 0`（如时变边界条件尚未激活）：$\Delta t_{\text{init}} \approx 0.01 \cdot \Delta t_{\max}$，保守但不过小。
- 当 `rate` 较大（系统快速演化）：$\Delta t_{\text{init}} \approx 0.01 \cdot d_0 / \text{rate}$，步长足够小以保证精度。

### 2.2 冷启动阶段 (`history_depth = 0`)

第一步无历史数据，无法构造预测器，因此：
- 仅执行 BE 求解，**不进行误差估计**。
- 下一步步长按最大增长因子放大：$\Delta t_{\text{next}} = \min(\Delta t \cdot \eta_{\max},\; \Delta t_{\max})$。

---

## 3. 预测器构造

### 3.1 线性预测器（1 阶，`history_depth ≥ 1`）

利用 $T^n$ 和 $T^{n-1}$ 两点做线性外推：

$$T_{\text{pred}}^{(1)} = (1 + r) \, T^n - r \, T^{n-1}, \quad r = \frac{\Delta t}{\Delta t_{\text{prev}}}$$

这是一阶 Taylor 展开的离散近似，截断误差为 $O(\Delta t^2)$。

### 3.2 二次 Lagrange 预测器（2 阶，`history_depth ≥ 2`）

利用 $T^n$、$T^{n-1}$、$T^{n-2}$ 三点构造二次 Lagrange 插值并外推到 $t^{n+1}$：

$$T_{\text{pred}}^{(2)} = L_0 \, T^{n-2} + L_1 \, T^{n-1} + L_2 \, T^n$$

其中 Lagrange 基函数值为：

$$L_0 = \frac{\Delta t \cdot (\Delta t + h_1)}{h_2 \cdot (h_1 + h_2)}, \quad
L_1 = \frac{-\Delta t \cdot (\Delta t + h_1 + h_2)}{h_1 \cdot h_2}, \quad
L_2 = \frac{(\Delta t + h_1) \cdot (\Delta t + h_1 + h_2)}{h_1 \cdot (h_1 + h_2)}$$

其中 $h_1 = \Delta t_{\text{prev}}$，$h_2 = \Delta t_{\text{prev2}}$。截断误差为 $O(\Delta t^3)$。

---

## 4. 校正器

### 4.1 后向欧拉（BE，1 阶校正器）

标准隐式欧拉离散：

$$(K + \frac{1}{\Delta t} C) \, T^{n+1} = F + \frac{1}{\Delta t} C \, T^n$$

局部截断误差为 $O(\Delta t^2)$。BE 解在每一步都会计算（作为 Picard 迭代的目标）。

### 4.2 变步长 BDF2（2 阶校正器）

对于不等间距时间步，BDF2 公式为：

$$\alpha_C \, T^{n+1} - \beta_n \, T^n + \beta_{n-1} \, T^{n-1} = C^{-1} (F - K \, T^{n+1})$$

定义 $r = \Delta t / \Delta t_{\text{prev}}$：

$$\alpha_C = \frac{1 + 2r}{(1+r)\,\Delta t}, \quad
\beta_n = \frac{1+r}{\Delta t}, \quad
\beta_{n-1} = \frac{r^2}{(1+r)\,\Delta t}$$

组装的线性系统：

$$(K + \alpha_C \, C) \, T^{n+1} = F + C \, (\beta_n \, T^n - \beta_{n-1} \, T^{n-1})$$

局部截断误差为 $O(\Delta t^3)$。

> **实现注意**：BDF2 的求解复用 BE 步中已组装的右端载荷 $F$（`F_current_`），因此 `SolveBDF2()` 必须在同一步的 `SolveBE()` 之后调用。系统矩阵 $(K + \alpha_C C)$ 按 `alpha_c` 值缓存，避免重复组装。

---

## 5. WRMS 误差范数

预测器与校正器之差的加权均方根范数：

$$\text{wrms} = \sqrt{\frac{1}{N_{\text{global}}} \sum_{i=1}^{N_{\text{global}}} \left(\frac{|T_{\text{corr},i} - T_{\text{pred},i}|}{\text{abstol} + \text{reltol} \cdot |T_{\text{corr},i}|}\right)^2}$$

**关键细节**：
- 求和使用 **真自由度**（True DOFs）：通过 `MPI_Allreduce` 对局部平方和求全局和，对局部自由度数求全局总数，避免共享节点被重复计算。
- 容差混合：每个自由度用 `abstol + reltol * |T_corr|` 加权，同时控制绝对误差和相对误差。
- **接受准则**：$\text{wrms} \leq 1.0$ 时接受该步。

---

## 6. 步长调整公式

$$\Delta t_{\text{new}} = \Delta t \cdot \eta_{\text{safety}} \cdot \text{wrms}^{-1/(p+1)}$$

其中 $p$ 为校正器阶数（1 或 2）。再施加限幅：

$$\Delta t_{\text{new}} = \text{clamp}(\Delta t_{\text{new}},\; \eta_{\min} \cdot \Delta t,\; \eta_{\max} \cdot \Delta t)$$

$$\Delta t_{\text{new}} = \text{clamp}(\Delta t_{\text{new}},\; \Delta t_{\min},\; \Delta t_{\max})$$

| 参数         | 含义                                       | 默认值 |
| ------------ | ------------------------------------------ | ------ |
| `eta_safety` | 安全因子，防止步长过于激进                 | 0.9    |
| `eta_max`    | 单步最大增长倍率                           | 2.0    |
| `eta_min`    | 单步最小收缩倍率                           | 0.2    |
| `dt_min`     | 绝对最小步长                               | 1e-6   |
| `dt_max`     | 绝对最大步长（0 = 使用 `output_interval`） | 0.0    |

**理论依据**：对于 $p$ 阶方法，局部截断误差 $\propto \Delta t^{p+1}$。若当前误差为 wrms，则要使误差降至 1 需步长缩减因子约 $\text{wrms}^{-1/(p+1)}$。因此：
- 1 阶：$\Delta t_{\text{new}} \propto \text{wrms}^{-1/2}$
- 2 阶：$\Delta t_{\text{new}} \propto \text{wrms}^{-1/3}$

---

## 7. 自适应阶数选择

当 `history_depth ≥ 2` 时，两种阶数同时计算各自的 wrms 和 $\Delta t_{\text{new}}$，按以下优先级选择：

```
if 两者都接受 (wrms ≤ 1):
    选择 dt_new 更大的阶数（效率优先）
else if 仅 2 阶接受:
    选 2 阶
else if 仅 1 阶接受:
    选 1 阶
else 两者都拒绝:
    取两者中更大的 dt_new 回退，拒绝该步
```

被选中阶数的校正器解将作为最终解。若选择 2 阶，则用 BDF2 解覆盖 BE 解：

```cpp
t_solver->GetTemperature().Distribute(t_solver->GetBDF2Solution());
```

---

## 8. 步拒绝与力接受机制

### 步拒绝

当 wrms > 1（两种阶数均不满足容差）时，拒绝当前步：
- 恢复温度场 `T = T_old`
- 用计算出的 `dt_new`（较小值）重试
- 累计 `rejected_steps` 计数

Picard 迭代不收敛时也触发步拒绝（步长减半重试）。

### 力接受 (`dt_min` 兜底)

当步长已降至 `dt_min`（判据：$\Delta t \leq \Delta t_{\min} \cdot (1 + 10^{-10})$）但仍不满足容差时，**强制接受**并输出警告日志：

```
Step force-accepted at dt_min=X.XXXXe-XX: wrms1=X.XXXXe+XX wrms2=X.XXXXe+XX
```

这保证算法不会无限缩小步长而陷入死循环。

---

## 9. 历史状态管理

### 状态变量

| 变量            | 含义                           |
| --------------- | ------------------------------ |
| `T_nm1_true`    | $T^{n-1}$ 的真自由度向量       |
| `T_nm2_true`    | $T^{n-2}$ 的真自由度向量       |
| `dt_prev`       | 上一步的 $\Delta t$            |
| `dt_prev2`      | 上上步的 $\Delta t$            |
| `history_depth` | 已保存的历史步数（0, 1, 或 2） |

### 更新逻辑

每步接受后，历史窗口滑动：

```cpp
T_nm2_true = T_nm1_true;       // n-2 ← n-1
T_nm1_true = T_n_true;          // n-1 ← n (当前步的 T_old)
dt_prev2 = dt_prev;
dt_prev = dt;
history_depth = min(history_depth + 1, 2);
```

`T_n_true`（当前步开始时的温度场真自由度）在误差估计阶段已计算，历史保存阶段直接复用，避免重复 `GetTrueDofs()` 调用。冷启动步（`history_depth=0`）无误差估计，则在历史保存阶段单独计算。

### 阶段转换

| `history_depth` |  可用预测器  | 可用校正器 | 说明               |
| :-------------: | :----------: | :--------: | ------------------ |
|        0        |      无      |     BE     | 冷启动，仅 BE 求解 |
|        1        | 线性（1 阶） |     BE     | 可估计 1 阶误差    |
|       2+        | 线性 + 二次  | BE + BDF2  | 完全自适应阶数选择 |

---

## 10. 输出快照插值

### 策略

采用 **线性插值** 在相邻时间步解之间生成输出快照。每完成一个时间步后，检查是否有输出时刻 $t_{\text{out}}$ 落入 $[t^n, t^{n+1}]$ 区间：

$$T(t_{\text{out}}) = (1 - \alpha) \, T(t^n) + \alpha \, T(t^{n+1}), \quad \alpha = \frac{t_{\text{out}} - t^n}{t^{n+1} - t^n}$$

输出插值在 **局部自由度**（Local DOFs）空间进行，因为 `ParGridFunction` 赋值需要匹配 `VSize()`。

> **设计选择**：经实测对比，线性插值在各测试用例上鲁棒性优于二次插值（二次插值在大步长比时可能引入振荡）。代码中保留了二次 Lagrange 插值函数 `EmitOutputSnapshot(3 snapshots)` 作为备用。

---

## 11. 时间循环总体流程

```
初始化: CacheMatrices(), 初始求解, ComputeInitialRate() → dt_init

while 尚有未输出的时刻:
    ① clamp(dt, dt_min, dt_max)
    ② 保存 T_old, EnableTransient(dt), 更新边界条件
    ③ Picard 迭代求解 BE
       - 不收敛 → 步拒绝（dt *= 0.5）或力接受
    ④ 预测-校正误差估计（跳过 history_depth=0）
       - 构造预测器（1 阶 / 2 阶）
       - 计算 wrms 和 dt_new
       - 自适应阶数选择
       - 不通过 → 步拒绝或力接受
       - 2 阶被选中 → BDF2 解覆盖 BE 解
    ⑤ 保存历史（T_nm1, T_nm2, dt_prev, dt_prev2）
    ⑥ t += dt, 保存当前快照
    ⑦ 输出落入当前区间的快照（线性插值）
    ⑧ dt = dt_new
```

---

## 12. 配置参数参考

```json
{
  "adaptive_dt": true,
  "adaptive_reltol": 1e-3,
  "adaptive_abstol": 1e-6,
  "dt_min": 1e-6,
  "dt_max": 0.0,
  "eta_safety": 0.9,
  "eta_max": 2.0,
  "eta_min": 0.2
}
```

### 参数调优建议

- **`adaptive_reltol`**：控制相对精度。典型范围 $10^{-2} \sim 10^{-5}$。过紧的容差并不一定改善输出精度（因输出依赖插值），但会增加步数。
- **`adaptive_abstol`**：控制绝对精度，防止小量场附近的分母趋零。应设为物理量的"噪声水平"。
- **`eta_safety`**：推荐 0.75–0.9。越小越保守（减少拒绝率，增加步数）。
- **`eta_max`**：推荐 2–5。限制步长突增，防止一步放大太多导致下一步拒绝。
- **`dt_min` / `dt_max`**：根据物理问题的时间尺度设置。`dt_max=0` 自动使用输出间隔。

---

## 13. 相关源文件

| 文件                                         | 职责                                               |
| -------------------------------------------- | -------------------------------------------------- |
| `src/transient/TransientSolver.cpp`          | 主时间循环、阶数选择、历史管理                     |
| `src/transient/TransientUtils.cpp`           | `ComputeAdaptiveError()`、`EmitOutputSnapshot()`   |
| `src/physics/ThermalFieldSolver.cpp`         | `SolveBE()`、`SolveBDF2()`、`ComputeInitialRate()` |
| `include/fem/transient/TransientSolver.hpp`  | `TransientSnapshot` 结构体、`TransientSolver` 类   |
| `include/fem/transient/TransientUtils.hpp`   | `WRMSResult` 结构体、工具函数声明                  |
| `include/fem/physics/ThermalFieldSolver.hpp` | `SolveState` 缓存结构、接口声明                    |

# design4 — 预测-校正自适应时间步长

本文档记录瞬态求解器从"简单自适应"升级为"预测-校正 + WRMS 误差控制"的完整设计。

---

## 1. 算法概述

原实现根据 Picard 迭代次数粗略调整步长（快收敛则放大 1.5×，慢收敛则缩小 0.5×），无法精确控制截断误差。新实现采用预测-校正（Predictor-Corrector）方案，通过 WRMS 误差范数在每步定量估计截断误差并据此控制步长。

完整流程为五步循环：

```
预测 → 校正(Picard) → 误差估计(WRMS) → 步长决策(accept/reject) → 前进或重试
```

---

## 2. 五步算法细节

### Step 1: 线性外推预测

利用前两步已知温度 $T_n, T_{n-1}$ 做一阶线性外推：

$$T_{n+1}^p = T_n + \Delta t \cdot \frac{T_n - T_{n-1}}{\Delta t_{\text{prev}}}$$

第一步（$n=0$）尚无历史，跳过预测，直接以初始 `transient_dt` 进行。

### Step 2: Picard 校正

以隐式 Backward-Euler 格式进行求解，内层 Picard 迭代耦合电热场直至收敛，得到修正值 $T_{n+1}$。

**Picard 容差要求**：必须比自适应步长容差小 1-2 个数量级。若 `adaptive_reltol = 1e-3`，则 `picard_tolerance` 应设为 `1e-5 ~ 1e-6`，否则 Picard 残差会污染截断误差估计，导致步长被错误缩小。

### Step 3: 误差估计

逐节点计算预测与修正的偏差 $e_i = |T_{n+1,i} - T_{n+1,i}^p|$。

### Step 4: WRMS 范数

$$\text{err} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(\frac{e_i}{\text{AbsTol} + \text{RelTol} \cdot |T_{n+1,i}|}\right)^2}$$

- **AbsTol**（绝对容差）：防止解接近零时分母失效。
- **RelTol**（相对容差）：控制精度要求，典型值 `1e-3`。
- $N$ 为**全局**自由度总数（MPI 并行时须 `MPI_Allreduce` 对 $\sum$ 和 $N$ 归约）。

设计约定：$\text{err} \le 1$ 表示步长可接受。

### Step 5: 步长控制

对于 Backward-Euler（$k=1$），最优步长公式：

$$\Delta t_{\text{new}} = \Delta t_{\text{old}} \cdot \text{safety} \cdot \left(\frac{1}{\text{err}}\right)^{1/(k+1)}$$

具体策略：

| 条件               | 动作                                                                                                                       |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------- |
| $\text{err} \le 1$ | **接受**，$\Delta t_{\text{new}} = \Delta t \cdot 0.9 \cdot (1/\text{err})^{1/2}$，增长因子限制在 $[0.5, 2.0]$             |
| $\text{err} > 1$   | **拒绝**，恢复 $T$ 到步前状态，$\Delta t_{\text{new}} = \Delta t \cdot 0.9 \cdot (1/\text{err})^{1/2}$，缩小因子不低于 0.2 |
| Picard 不收敛      | **拒绝**，$\Delta t \times 0.5$，重试                                                                                      |

步长始终被裁剪到 $[\Delta t_{\min}, \Delta t_{\max}]$，其中 $\Delta t_{\max}$ 默认等于 `transient_output_interval`。

---

## 3. 输出插值

自适应步长下实际求解时刻不落在均匀输出网格上。因此：

1. **求解阶段**：每个被接受的时间步都保存原始快照（raw snapshot）。
2. **导出阶段**：对原始快照序列在均匀 `transient_output_interval` 网格上做**线性插值**，生成均匀快照后再导出。

插值公式（在 $[t_j, t_{j+1}]$ 区间内对任一场 $\mathbf{u}$）：

$$\mathbf{u}(t_{\text{out}}) = (1-\alpha)\,\mathbf{u}(t_j) + \alpha\,\mathbf{u}(t_{j+1}), \quad \alpha = \frac{t_{\text{out}} - t_j}{t_{j+1} - t_j}$$

实现为 `MultiPhysicsCoupler::InterpolateSnapshots()`，对 voltage / temperature / displacement 三个场向量分别插值。

---

## 4. MPI 并行一致性

在多进程并行下，每个 rank 仅持有本地 DOF 子集。两个关键决策点必须进行全局归约，否则各 rank 做出不同的 accept/reject 判断，导致 MPI 通信不匹配而崩溃（`MPI_ERR_TRUNCATE`）：

### Picard 收敛检查

```cpp
double diff_sq = delta * delta;   // 本地 dot product
double ref_sq  = T_prev * T_prev;
#ifdef MFEM_USE_MPI
if (config_.fe.IsParallel())
{
    double buf[2] = {diff_sq, ref_sq}, gbuf[2];
    MPI_Allreduce(buf, gbuf, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    diff_sq = gbuf[0];  ref_sq = gbuf[1];
}
#endif
double rel_change = std::sqrt(diff_sq) / std::sqrt(ref_sq);
```

### WRMS 误差

```cpp
double local_sum_sq = ...;  // 本地加权平方和
int    local_N      = ...;  // 本地 DOF 数
double global_sum_sq;  long long global_N;
#ifdef MFEM_USE_MPI
if (config_.fe.IsParallel())
{
    MPI_Allreduce(&local_sum_sq, &global_sum_sq, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_N_ll,  &global_N,      1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
}
#endif
wrms_err = std::sqrt(global_sum_sq / global_N);
```

所有 rank 使用相同的 `wrms_err` 值做 accept/reject 决策，保证代码路径一致。

---

## 5. 配置接口

### 新增 JSON 字段（`simulation` 节）

| 字段              | 类型   | 默认值 | 说明                                                 |
| ----------------- | ------ | ------ | ---------------------------------------------------- |
| `adaptive_dt`     | bool   | false  | 启用预测-校正自适应步长                              |
| `adaptive_reltol` | double | 1e-3   | WRMS 相对容差                                        |
| `adaptive_abstol` | double | 1e-6   | WRMS 绝对容差                                        |
| `dt_min`          | double | 1e-6   | 最小允许步长                                         |
| `dt_max`          | double | 0      | 最大允许步长（0 = 使用 `transient_output_interval`） |

### Config.hpp 对应结构

```cpp
struct SimulationConfig {
    // ... existing fields ...
    bool   adaptive_dt      = false;
    double adaptive_reltol  = 1.0e-3;
    double adaptive_abstol  = 1.0e-6;
    double dt_min           = 1.0e-6;
    double dt_max           = 0.0;
};
```

### 推荐配置示例

```json
{
  "picard_tolerance": 1e-6,
  "adaptive_dt": true,
  "adaptive_reltol": 1e-3,
  "adaptive_abstol": 1e-6,
  "dt_min": 1e-3,
  "dt_max": 20.0
}
```

注意 `picard_tolerance`（1e-6）比 `adaptive_reltol`（1e-3）小 3 个数量级，满足精度隔离要求。

---

## 6. 代码改动清单

| 文件                                           | 改动                                                                                                                      |
| ---------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `include/fem/frontend/Config.hpp`              | `SimulationConfig` 增加 5 个自适应字段                                                                                    |
| `src/frontend/ConfigLoader.cpp`                | 解析 `adaptive_reltol / adaptive_abstol / dt_min / dt_max`                                                                |
| `include/fem/coupling/MultiPhysicsCoupler.hpp` | 新增 `InterpolateSnapshots()` 声明                                                                                        |
| `src/coupling/MultiPhysicsCoupler.cpp`         | **重写 `SolveTransient()`**（预测-校正-WRMS 循环）；新增 `InterpolateSnapshots()` 实现；Picard 和 WRMS 加 `MPI_Allreduce` |
| `docs/physics.schema.json`                     | 新增 4 个字段定义，更新 `adaptive_dt` 描述                                                                                |
| `configs/busbar_etm_nonlinear_transient.json`  | 添加自适应参数                                                                                                            |

---

## 7. 验证结果

Busbar 电热力非线性瞬态（t=0→100s），2 进程 / 4 进程均通过：

```
=== Transient solve: t=[0, 100], dt=10, output_interval=10 ===
  Adaptive dt: reltol=1.00e-03, abstol=1.00e-06, dt_range=[1.00e-03, 2.00e+01]

Step 1: t=10.0000, dt=1.0000e+01, Picard=2        (首步无预测)
Step rejected: WRMS err=1.16e+00 > 1, new dt=8.35e+00
Step 2: t=18.3544, dt=8.35e+00, Picard=2, WRMS=8.97e-01
Step 3: t=26.2956, dt=7.94e+00, Picard=2, WRMS=3.84e-01
Step 4: t=37.8330, dt=1.15e+01, Picard=2, WRMS=3.18e-01
Step 5: t=56.2476, dt=1.84e+01, Picard=2, WRMS=2.40e-01
Step 6: t=76.2476, dt=2.00e+01, Picard=2, WRMS=8.62e-02  (触顶 dt_max)
Step 7: t=96.2476, dt=2.00e+01, Picard=2, WRMS=3.42e-02
Step 8: t=100.000, dt=3.75e+00, Picard=2, WRMS=1.01e-03  (末步精确到达 t_end)

8 accepted steps, 1 rejected, 9 raw snapshots → 11 interpolated output snapshots
```

行为符合预期：

- 首步用固定 `transient_dt=10`，第二步起预测-校正生效。
- WRMS > 1 时正确拒绝并缩小步长。
- 随温度趋于稳态，WRMS 降低，步长自动增长直至触顶 `dt_max=20`。
- 最终 11 个均匀快照（t=0,10,...,100）通过线性插值从 9 个不均匀原始快照生成。

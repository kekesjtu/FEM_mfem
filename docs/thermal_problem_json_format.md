# 我们项目的热问题 JSON 格式

本格式用于在当前 FEM_mfem 项目中描述稳态热传导问题，支持：
- 域分段材料导热系数
- 域分段热源
- Dirichlet 与 Robin 边界

## 顶层结构

```json
{
  "simulation": { ... },
  "materials": [ ... ],
  "domain_materials": { ... },
  "fields": [ ... ]
}
```

## simulation

- `mesh_path`：网格文件路径（当前要求 MFEM 可读，如 `.msh`）
- `order`：有限元阶次
- `uniform_refinement_levels`：一致加密层数
- `output_dir`：输出目录
- `log_level`：日志级别

## materials

每个材料定义属性表达式，热问题至少建议给出 `diffusion`（导热系数）。

```json
{ "name": "Aluminum", "properties": { "diffusion": "237.0" } }
```

## domain_materials

将材料映射到域属性号（mesh volume attributes）：

```json
"domain_materials": {
  "Aluminum": [3,4,5],
  "Silicon": [11]
}
```

程序会自动按该映射构造分段 `diffusion` 系数。

## fields[*]

热场字段推荐配置：

- `name`：场名（如 `temperature`）
- `type`：`thermal`
- `source_default`：默认体热源
- `source_by_domain`：按域属性覆盖热源值
- `dirichlet_value`：Dirichlet 温度值
- `dirichlet_bdr_attributes`：Dirichlet 边界属性号
- `robin_l`：Robin 左端系数 $L$（方程形式：$k\partial_n T + L T = q$）
- `robin_q`：Robin 右端项 $q$
- `robin_on_remaining_boundaries`：是否将 Robin 施加到“除 Dirichlet 外的所有边界”

## 示例

完整示例见 [configs/complex_thermal_robin_project.json](../configs/complex_thermal_robin_project.json)。

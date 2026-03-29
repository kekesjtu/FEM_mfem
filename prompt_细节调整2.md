## 一、机械场边界与载荷能力补全  
1. 目前项目不仅需要支持分域的 normal_displacement_boundary_conditions。  
2. 同时也需要支持分域的 displacement_boundary_conditions。  
3. 还需要支持分域体力 body_force_by_domain。  
4. body_force_default 作为默认值，处理方式应与 source_default 一致（先设默认，再按分域覆盖）。

## 二、配置层自动补全拉梅参数条目  
1. 在 config 中增加自动补全条目的功能。  
2. 当识别到 young_modulus 与 poisson_ratio 时，自动向 material database 增加拉梅系数相关条目。  
3. 完成后可删除下面这段“从杨氏模量和泊松比转换拉梅参数”的逻辑。  
4. 后续统一使用项目内的 PiecewiseConstantCoefficient。

```cpp
可删除的现有转换逻辑：
    // Build Lame parameters from Young's modulus and Poisson ratio
    int max_attr = coeff::GetMaxAttribute(mesh);
    mfem::PWConstCoefficient lambda_pw(max_attr);
    mfem::PWConstCoefficient mu_pw(max_attr);

    std::unordered_map<int, std::string> attr_to_mat;
    for (const auto &[mat_name, domains] : config_.materials.domain_to_material)
    {
        for (int a : domains)
        {
            attr_to_mat[a] = mat_name;
        }
    }

    for (int a = 1; a <= max_attr; ++a)
    {
        auto it = attr_to_mat.find(a);
        if (it == attr_to_mat.end())
        {
            lambda_pw(a) = 0.0;
            mu_pw(a) = 0.0;
            continue;
        }

        const auto &props = config_.materials.material_to_properties.at(it->second);

        auto E_it = props.find("young_modulus");
        auto nu_it = props.find("poisson_ratio");
        if (E_it == props.end() || nu_it == props.end())
        {
            lambda_pw(a) = 0.0;
            mu_pw(a) = 0.0;
            continue;
        }

        double E = std::stod(E_it->second);
        double nu = std::stod(nu_it->second);
        lambda_pw(a) = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        mu_pw(a) = E / (2.0 * (1.0 + nu));

        logger->debug("Attr {} ({}): E={:.3e}, nu={:.4f}, lambda={:.3e}, mu={:.3e}", a, it->second,
                      E, nu, lambda_pw(a), mu_pw(a));
    }
```

## 三、并行能力建设与验证  
1. 加入完整并行功能，参考 mfem 并行案例。  
2. 详细梳理并正确使用 mfem 并行计算流程。  
3. 验证并行速度与收敛性是否符合预期。  
4. 采用条件编译区分并行与串行代码路径。  
5. 串行与并行两套实现保持同一套对外接口。

## 四、优化脚本
优化脚本compare_solution_txt，更好地适配位移向量场的比较。

以上内容完成后,交付文档design3
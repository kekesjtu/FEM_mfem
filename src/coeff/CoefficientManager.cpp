#include "coeff/CoefficientManager.hpp"

#include "log/Logger.hpp"

#include <atomic>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace fem::coeff
{
namespace
{
std::atomic_size_t next_expression_coefficient_runtime_cache_id{1};

std::size_t AllocateExpressionCoefficientRuntimeCacheId()
{
    return next_expression_coefficient_runtime_cache_id.fetch_add(1, std::memory_order_relaxed);
}

struct
    ExpressionCoefficientRuntimeCache  // 这个结构体会在每个 ExpressionCoefficient 的 Eval 中以
                                       // thread_local 方式使用，缓存该 ExpressionCoefficient
                                       // 在当前线程中的运行时数据，包括编译后的表达式和物理坐标（不用每次new）
{
    std::vector<frontend::ExpressionRuntime> expr_runtimes;
    mfem::Vector phys_point;
};
}  // namespace

bool IsPropertyConstant(const std::string &property_name, const frontend::MaterialDatabase &db)
{
    for (const auto &[attr, props] : db.attr_to_properties)
    {
        auto it = props.find(property_name);
        if (it != props.end())
        {
            // Check if it's a pure number
            const std::string &expr_str = it->second;
            try
            {
                size_t processed = 0;
                std::stod(expr_str, &processed);
                if (processed != expr_str.size())
                {
                    return false;
                }
            }
            catch (...)
            {
                return false;
            }
        }
    }
    return true;
}

// --- ExpressionCoefficient ---

ExpressionCoefficient::ExpressionCoefficient(const std::string &property_name,
                                             const frontend::MaterialDatabase &db, mfem::Mesh &mesh,
                                             mfem::GridFunction *const *temperature_gf)
    : temperature_gf_(temperature_gf),
      mesh_(mesh),
      runtime_cache_id_(
          AllocateExpressionCoefficientRuntimeCacheId())  // 每个 ExpressionCoefficient
                                                          // 构造时领一个唯一 ID（1, 2, 3...），存入
                                                          // runtime_cache_id_
{
    const int max_attr = mesh.attributes.Max();
    // Default expression: "0.0" for unassigned attributes.
    attr_to_expr_.resize(max_attr + 1,
                         AttributeExpression{frontend::Expression("0.0"), true, 0.0});

    // For each attribute id, look up the property expression prepared by LoadMaterials().
    for (int attr = 1; attr <= max_attr; ++attr)
    {
        auto prop_it = db.attr_to_properties.find(attr);
        if (prop_it == db.attr_to_properties.end())
        {
            continue;
        }

        auto expr_it = prop_it->second.find(property_name);
        if (expr_it == prop_it->second.end())
        {
            continue;
        }

        frontend::Expression expr(expr_it->second);
        // Check if constant
        bool is_const = false;
        double const_val = 0.0;
        try
        {
            size_t processed = 0;
            const_val = std::stod(expr_it->second, &processed);
            is_const = (processed == expr_it->second.size());
        }
        catch (...)
        {
            is_const = false;
        }

        attr_to_expr_[attr] = AttributeExpression{std::move(expr), is_const, const_val};
    }
}

double ExpressionCoefficient::Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
{
    const int attr = T.Attribute;
    if (attr < 0 || attr >= static_cast<int>(attr_to_expr_.size()))
    {
        return 0.0;
    }

    auto &de = attr_to_expr_[attr];
    if (de.is_constant)
    {
        return de.constant_value;
    }

    // Need to evaluate expression with spatial coords and temperature
    thread_local std::unordered_map<std::size_t, ExpressionCoefficientRuntimeCache> runtime_caches;
    auto &cache = runtime_caches[runtime_cache_id_];
    if (cache.expr_runtimes.size() != attr_to_expr_.size())
    {
        cache.expr_runtimes.clear();
        cache.expr_runtimes.reserve(attr_to_expr_.size());
        for (const auto &entry : attr_to_expr_)
        {
            cache.expr_runtimes.emplace_back(entry.expr);
        }
    }

    frontend::EvalContext ctx;
    cache.phys_point.SetSize(T.GetSpaceDim());
    T.Transform(ip, cache.phys_point);
    ctx.x = cache.phys_point.Size() > 0 ? cache.phys_point(0) : 0.0;
    ctx.y = cache.phys_point.Size() > 1 ? cache.phys_point(1) : 0.0;
    ctx.z = cache.phys_point.Size() > 2 ? cache.phys_point(2) : 0.0;

    if (temperature_gf_ && *temperature_gf_)
    {
        ctx.T = (*temperature_gf_)->GetValue(T, ip);
    }
    else
    {
        ctx.T = reference_temperature_;
    }

    return cache.expr_runtimes[attr].Evaluate(ctx);
}

// --- PiecewiseConstantCoefficient ---

PiecewiseConstantCoefficient::PiecewiseConstantCoefficient(const std::string &property_name,
                                                           const frontend::MaterialDatabase &db,
                                                           mfem::Mesh &mesh)
{
    const int max_attr = mesh.attributes.Max();
    attr_to_value_.resize(max_attr + 1, 0.0);

    auto logger = fem::log::Get();
    for (int attr = 1; attr <= max_attr; ++attr)
    {
        auto prop_it = db.attr_to_properties.find(attr);
        if (prop_it == db.attr_to_properties.end())
        {
            continue;
        }

        auto expr_it = prop_it->second.find(property_name);
        if (expr_it == prop_it->second.end())
        {
            continue;
        }

        try
        {
            attr_to_value_[attr] = std::stod(expr_it->second);
        }
        catch (const std::exception &e)
        {
            std::string message = "PiecewiseConstantCoefficient: property '";
            message += property_name;
            message += "' for material '";
            auto material_it = db.attr_to_material.find(attr);
            message += material_it != db.attr_to_material.end() ? material_it->second : "<unknown>";
            message += "' is not a constant: ";
            message += expr_it->second;
            throw std::runtime_error(message);
        }

        auto material_it = db.attr_to_material.find(attr);
        const std::string material_name =
            material_it != db.attr_to_material.end() ? material_it->second : "<unknown>";
        logger->debug("Material '{}' attr={}: {}={}", material_name, attr, property_name,
                      attr_to_value_[attr]);
    }
}

double PiecewiseConstantCoefficient::Eval(mfem::ElementTransformation &T,
                                          const mfem::IntegrationPoint & /*ip*/)
{
    const int attr = T.Attribute;
    if (attr < 0 || attr >= static_cast<int>(attr_to_value_.size()))
    {
        return 0.0;
    }
    return attr_to_value_[attr];
}

// --- JouleHeatingCoefficient ---

JouleHeatingCoefficient::JouleHeatingCoefficient(mfem::Coefficient *const *sigma,
                                                 mfem::GridFunction *const *voltage)
    : sigma_(sigma), voltage_(voltage)
{
}

double JouleHeatingCoefficient::Eval(mfem::ElementTransformation &T,
                                     const mfem::IntegrationPoint &ip)
{
    if (!sigma_ || !*sigma_ || !voltage_ || !*voltage_)
    {
        throw std::runtime_error(
            "JouleHeatingCoefficient requires both conductivity and voltage fields before Eval");
    }

    // Q = sigma * |grad(V)|^2
    const double sigma_val = (*sigma_)->Eval(T, ip);

    const int dim = T.GetSpaceDim();
    thread_local mfem::Vector grad_v_buf;
    grad_v_buf.SetSize(dim);
    (*voltage_)->GetGradient(T, grad_v_buf);

    double grad_sq = 0.0;
    for (int d = 0; d < dim; ++d)
    {
        grad_sq += grad_v_buf(d) * grad_v_buf(d);
    }

    return sigma_val * grad_sq;
}

double CombinedHeatSourceCoefficient::Eval(mfem::ElementTransformation &T,
                                           const mfem::IntegrationPoint &ip)
{
    double value = 0.0;

    if (base_source_)
        value += base_source_->Eval(T, ip);

    if (joule_source_)
        value += joule_source_->Eval(T, ip);

    return value;
}

}  // namespace fem::coeff

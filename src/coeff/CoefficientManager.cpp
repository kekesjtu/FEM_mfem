#include "fem/coeff/Coeffmanagaer.hpp"

#include "fem/log/Logger.hpp"

#include <cmath>
#include <stdexcept>

namespace fem::coeff
{

int GetMaxAttribute(mfem::Mesh &mesh)
{
    int max_attr = 0;
    for (int i = 0; i < mesh.GetNE(); ++i)
    {
        max_attr = std::max(max_attr, mesh.GetAttribute(i));
    }
    return max_attr;
}

bool IsPropertyConstant(const std::string &property_name, const frontend::MaterialDatabase &db)
{
    for (const auto &[mat_name, props] : db.material_to_properties)
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
                                             mfem::GridFunction *temperature_gf)
    : temperature_gf_(temperature_gf), mesh_(mesh)
{
    const int max_attr = GetMaxAttribute(mesh);
    // Default expression: "0.0" for unassigned domains
    attr_to_expr_.resize(max_attr + 1, DomainExpression{frontend::Expression("0.0"), true, 0.0});

    // Build attribute -> material name map
    std::unordered_map<int, std::string> attr_to_mat;
    for (const auto &[mat_name, domains] : db.domain_to_material)
    {
        for (int attr : domains)
        {
            attr_to_mat[attr] = mat_name;
        }
    }

    // For each attribute, look up the property expression
    for (int attr = 1; attr <= max_attr; ++attr)
    {
        auto mat_it = attr_to_mat.find(attr);
        if (mat_it == attr_to_mat.end())
        {
            continue;
        }

        const std::string &mat_name = mat_it->second;
        auto prop_it = db.material_to_properties.find(mat_name);
        if (prop_it == db.material_to_properties.end())
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

        attr_to_expr_[attr] = DomainExpression{std::move(expr), is_const, const_val};
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
    frontend::EvalContext ctx;
    mfem::Vector phys_point;
    T.Transform(ip, phys_point);
    ctx.x = phys_point.Size() > 0 ? phys_point(0) : 0.0;
    ctx.y = phys_point.Size() > 1 ? phys_point(1) : 0.0;
    ctx.z = phys_point.Size() > 2 ? phys_point(2) : 0.0;

    if (temperature_gf_)
    {
        ctx.T = temperature_gf_->GetValue(T, ip);
    }
    else
    {
        ctx.T = reference_temperature_;
    }

    return de.expr.Evaluate(ctx);
}

// --- PiecewiseConstantCoefficient ---

PiecewiseConstantCoefficient::PiecewiseConstantCoefficient(const std::string &property_name,
                                                           const frontend::MaterialDatabase &db,
                                                           mfem::Mesh &mesh)
{
    const int max_attr = GetMaxAttribute(mesh);
    attr_to_value_.resize(max_attr + 1, 0.0);

    std::unordered_map<int, std::string> attr_to_mat;
    for (const auto &[mat_name, domains] : db.domain_to_material)
    {
        for (int attr : domains)
        {
            attr_to_mat[attr] = mat_name;
        }
    }

    auto logger = fem::log::Get();
    for (int attr = 1; attr <= max_attr; ++attr)
    {
        auto mat_it = attr_to_mat.find(attr);
        if (mat_it == attr_to_mat.end())
        {
            continue;
        }

        const std::string &mat_name = mat_it->second;
        auto prop_it = db.material_to_properties.find(mat_name);
        if (prop_it == db.material_to_properties.end())
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
            throw std::runtime_error("PiecewiseConstantCoefficient: property '" + property_name +
                                     "' for material '" + mat_name +
                                     "' is not a constant: " + expr_it->second);
        }

        logger->debug("Material '{}' attr={}: {}={}", mat_name, attr, property_name,
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

JouleHeatingCoefficient::JouleHeatingCoefficient(mfem::Coefficient &sigma,
                                                 mfem::GridFunction &voltage)
    : sigma_(sigma), voltage_(voltage)
{
}

double JouleHeatingCoefficient::Eval(mfem::ElementTransformation &T,
                                     const mfem::IntegrationPoint &ip)
{
    // Q = sigma * |grad(V)|^2
    const double sigma_val = sigma_.Eval(T, ip);

    const int dim = T.GetSpaceDim();
    mfem::Vector grad_v(dim);
    voltage_.GetGradient(T, grad_v);

    double grad_sq = 0.0;
    for (int d = 0; d < dim; ++d)
    {
        grad_sq += grad_v(d) * grad_v(d);
    }

    return sigma_val * grad_sq;
}

}  // namespace fem::coeff

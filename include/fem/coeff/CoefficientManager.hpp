// CoefficientManager: bridges Config material data to MFEM Coefficient objects.
// Supports piecewise-constant and piecewise-expression (temperature-dependent) coefficients.
#pragma once

#include "fem/frontend/Config.hpp"
#include "fem/frontend/Expression.hpp"
#include "mfem.hpp"

#include <string>
#include <vector>

namespace fem::coeff
{

/// GridFunctionCoefficient that evaluates a material expression at integration points,
/// optionally depending on a temperature GridFunction for coupled problems.
class ExpressionCoefficient : public mfem::Coefficient
{
  public:
    /// Build a piecewise expression coefficient from material database.
    /// @param property_name  e.g. "electrical_conductivity", "diffusion"
    /// @param db             material database with domain->material and material->properties maps
    /// @param mesh           the mesh (to get element attributes)
    /// @param temperature_gf optional late-bound temperature field handle for T-dependent
    /// expressions
    ExpressionCoefficient(const std::string &property_name, const frontend::MaterialDatabase &db,
                          mfem::Mesh &mesh, mfem::GridFunction * const *temperature_gf = nullptr);

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;

    /// Set the fallback temperature (K) used when no temperature GridFunction is available.
    void SetReferenceTemperature(double T_ref)
    {
        reference_temperature_ = T_ref;
    }

  private:
    struct DomainExpression
    {
        frontend::Expression expr;
        bool is_constant;
        double constant_value;
    };

    // attribute -> expression (1-based indexing, slot 0 unused)
    std::vector<DomainExpression> attr_to_expr_;
    mfem::GridFunction * const *temperature_gf_;
    mfem::Mesh &mesh_;
    double reference_temperature_ = 293.15;
    mutable mfem::Vector phys_point_buf_;  // pre-allocated buffer for Transform
};

/// Simple piecewise constant coefficient built from material property.
/// Faster path when no expression depends on T/x/y/z.
class PiecewiseConstantCoefficient : public mfem::Coefficient
{
  public:
    PiecewiseConstantCoefficient(const std::string &property_name,
                                 const frontend::MaterialDatabase &db, mfem::Mesh &mesh);

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;

  private:
    std::vector<double> attr_to_value_;  // 1-based
};

/// Joule heating source: Q = sigma * |grad(V)|^2, computed at integration points.
class JouleHeatingCoefficient : public mfem::Coefficient
{
  public:
    JouleHeatingCoefficient(mfem::Coefficient * const *sigma = nullptr,
                            mfem::GridFunction * const *voltage = nullptr);

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;

  private:
    mfem::Coefficient * const *sigma_ = nullptr;
    mfem::GridFunction * const *voltage_ = nullptr;
    mutable mfem::Vector grad_v_buf_;  // pre-allocated buffer for gradient
};

class CombinedHeatSourceCoefficient : public mfem::Coefficient
{
  public:
    CombinedHeatSourceCoefficient(mfem::Coefficient *base_source,
                                  JouleHeatingCoefficient *joule_source)
        : base_source_(base_source), joule_source_(joule_source)
    {
    }

    double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override;

  private:
    mfem::Coefficient *base_source_ = nullptr;
    JouleHeatingCoefficient *joule_source_ = nullptr;
};

/// Utility: check if a material property is purely constant across all domains
bool IsPropertyConstant(const std::string &property_name, const frontend::MaterialDatabase &db);

}  // namespace fem::coeff
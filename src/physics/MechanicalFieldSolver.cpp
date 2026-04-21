#include "fem/physics/MechanicalFieldSolver.hpp"

#include "fem/assembly/MechanicalAssembler.hpp"
#include "fem/coeff/CoefficientManager.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

namespace fem::physics
{

MechanicalFieldSolver::MechanicalFieldSolver(frontend::ProjectConfig &config)
    : config_(config), displacement_(&config.fe.GetVectorFESpace())
{
    displacement_ = 0.0;

    auto &mesh = config_.fe.GetMesh();

    // Cache constant coefficients (lambda, mu, body_force, alpha)
    cached_lambda_ = std::make_unique<coeff::PiecewiseConstantCoefficient>("lambda_lame",
                                                                           config_.materials, mesh);
    cached_mu_ =
        std::make_unique<coeff::PiecewiseConstantCoefficient>("mu_lame", config_.materials, mesh);

    const auto &field_config = config_.mechanical_field;
    const int dim = mesh.Dimension();

    if (!field_config.domain_to_body_force.empty())
    {
        int max_attr = coeff::GetMaxAttribute(mesh);
        auto pw_vec = std::make_unique<mfem::VectorArrayCoefficient>(dim);
        for (int d = 0; d < dim; ++d)
        {
            mfem::Vector comp_vals(max_attr);
            double default_val = (d < static_cast<int>(field_config.body_force_default.size()))
                                     ? field_config.body_force_default[d]
                                     : 0.0;
            for (int a = 0; a < max_attr; ++a)
            {
                auto it = field_config.domain_to_body_force.find(a + 1);
                if (it != field_config.domain_to_body_force.end() &&
                    d < static_cast<int>(it->second.size()))
                    comp_vals(a) = it->second[d];
                else
                    comp_vals(a) = default_val;
            }
            pw_vec->Set(d, new mfem::PWConstCoefficient(comp_vals));
        }
        cached_body_force_ = std::move(pw_vec);
    }
    else
    {
        mfem::Vector body_force_vec(dim);
        body_force_vec = 0.0;
        for (int d = 0; d < dim && d < static_cast<int>(field_config.body_force_default.size());
             ++d)
            body_force_vec(d) = field_config.body_force_default[d];
        cached_body_force_ = std::make_unique<mfem::VectorConstantCoefficient>(body_force_vec);
    }

    // alpha is cached but temp_coeff is per-solve (depends on current temperature GF)
    // alpha will be built when temperature_gf_ is set via SetTemperatureField

    // Build and cache BC markers once
    BuildBCMarkers();
}

void MechanicalFieldSolver::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    temperature_gf_ = temperature_gf;

    // Build alpha coefficient now that we know we have thermal coupling
    if (temperature_gf_ && !cached_alpha_)
    {
        cached_alpha_ = std::make_unique<coeff::PiecewiseConstantCoefficient>(
            "thermal_expansion_secant", config_.materials, config_.fe.GetMesh());
    }
}

void MechanicalFieldSolver::BuildBCMarkers()
{
    auto &mesh = config_.fe.GetMesh();
    const auto &field_config = config_.mechanical_field;
    const int num_bdr = mesh.bdr_attributes.Max();
    const int dim = mesh.Dimension();

    // Displacement (essential) BCs
    cached_essential_bdr_.SetSize(num_bdr);
    cached_essential_bdr_ = 0;
    cached_disp_coeffs_.reserve(field_config.displacement_bcs.size());
    cached_disp_markers_.reserve(field_config.displacement_bcs.size());
    for (const auto &dbc : field_config.displacement_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : dbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
            {
                marker[attr - 1] = 1;
                cached_essential_bdr_[attr - 1] = 1;
            }
        }
        mfem::Vector val(dim);
        val = 0.0;
        for (int d = 0; d < dim && d < static_cast<int>(dbc.value.size()); ++d)
            val(d) = dbc.value[d];
        cached_disp_coeffs_.emplace_back(val);
        cached_disp_markers_.push_back(std::move(marker));
    }

    // Normal displacement penalty BCs
    cached_penalty_coeffs_.reserve(field_config.normal_displacement_bcs.size());
    cached_nd_markers_.reserve(field_config.normal_displacement_bcs.size());
    for (const auto &nbc : field_config.normal_displacement_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : nbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        cached_penalty_coeffs_.emplace_back(nbc.penalty);
        cached_nd_markers_.push_back(std::move(marker));
    }

    // Pressure BCs
    cached_pressure_coeffs_.reserve(field_config.pressure_bcs.size());
    cached_pressure_markers_.reserve(field_config.pressure_bcs.size());
    for (const auto &pbc : field_config.pressure_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : pbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        cached_pressure_coeffs_.emplace_back(pbc.value);
        cached_pressure_markers_.push_back(std::move(marker));
    }

    // Cache essential tdofs
    config_.fe.GetVectorFESpace().GetEssentialTrueDofs(cached_essential_bdr_,
                                                       cached_essential_tdofs_);
}

std::unique_ptr<mfem::Coefficient> MechanicalFieldSolver::BuildTempCoefficient()
{
    if (temperature_gf_)
        return std::make_unique<mfem::GridFunctionCoefficient>(temperature_gf_);
    return nullptr;
}

void MechanicalFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->debug("Mechanical full solve");

    // Use CacheStiffnessMatrix + SolveRHSOnly pattern to avoid redundant assembly/factorization
    CacheStiffnessMatrix();
    SolveRHSOnly();
}

// ------------------------------------------------------------------
// CacheStiffnessMatrix — delegate to MechanicalAssembler
// ------------------------------------------------------------------
void MechanicalFieldSolver::CacheStiffnessMatrix()
{
    if (stiffness_cached_)
        return;

    auto logger = fem::log::Get();
    logger->debug("Caching mechanical stiffness matrix...");

    cached_bilinear_ = assembly::MechanicalAssembler::AssembleStiffness(
        config_.fe.GetVectorFESpace(), *cached_lambda_, *cached_mu_, cached_penalty_coeffs_,
        cached_nd_markers_);

    stiffness_cached_ = true;
    solver_factorized_ = false;
    cached_solver_.reset();
    logger->debug("Stiffness matrix cached.");
}

// ------------------------------------------------------------------
// SolveRHSOnly — delegate RHS assembly to MechanicalAssembler, solve with cached K
// ------------------------------------------------------------------
void MechanicalFieldSolver::SolveRHSOnly()
{
    auto logger = fem::log::Get();
    logger->debug("Mechanical solve (cached K, RHS-only)");

    if (!stiffness_cached_)
    {
        logger->warn("Stiffness matrix not cached, falling back to full Solve()");
        Solve();
        return;
    }

    // Only rebuild the temperature-dependent coefficient
    auto temp_coeff = BuildTempCoefficient();

    // Assemble RHS via assembler
    auto lf = assembly::MechanicalAssembler::AssembleRHS(
        config_.fe.GetVectorFESpace(), *cached_body_force_, cached_pressure_coeffs_,
        cached_pressure_markers_, cached_lambda_.get(), cached_mu_.get(), cached_alpha_.get(),
        temp_coeff.get(), config_.mechanical_field.reference_temperature);

    // Apply essential BCs to initial guess
    displacement_ = 0.0;
    for (size_t i = 0; i < cached_disp_markers_.size(); ++i)
        displacement_.ProjectBdrCoefficient(cached_disp_coeffs_[i], cached_disp_markers_[i]);

    // Form linear system using cached bilinear form
    mfem::OperatorPtr A;
    mfem::Vector X, B;
    cached_bilinear_->FormLinearSystem(cached_essential_tdofs_, displacement_, *lf, A, X, B);

    if (!solver_factorized_)
    {
        cached_solver_ = solver::CreateLinearSolver(config_.simulation.GetSolver("mechanical"),
                                                    1e-12, 1e-12, 5000, 0);
        cached_solver_->Solve(*A, B, X);
        solver_factorized_ = true;
    }
    else
    {
        cached_solver_->Mult(B, X);
    }
    cached_bilinear_->RecoverFEMSolution(X, *lf, displacement_);

    logger->debug("Mechanical: U=[{:.6e}, {:.6e}]", displacement_.Min(), displacement_.Max());
}

}  // namespace fem::physics

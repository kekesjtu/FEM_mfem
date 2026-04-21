#include "fem/physics/ElectrostaticFieldSolver.hpp"

#include "fem/coeff/CoefficientManager.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

namespace fem::physics
{

ElectrostaticFieldSolver::ElectrostaticFieldSolver(frontend::ProjectConfig &config)
    : config_(config), voltage_(&config.fe.GetScalarFESpace())
{
    voltage_ = 0.0;

    // Pre-parse time-varying BC expressions
    for (const auto &dbc : config_.electric_field.dirichlet_bcs)
    {
        if (!dbc.value_expr.empty())
            bc_expressions_.emplace_back(dbc.value_expr);
        else
            bc_expressions_.emplace_back(std::to_string(dbc.value));
    }

    // Build and cache BC markers once (they never change)
    BuildBCMarkers();
}

void ElectrostaticFieldSolver::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    temperature_gf_ = temperature_gf;
}

void ElectrostaticFieldSolver::SetElectricalConductivity(mfem::Coefficient *sigma)
{
    sigma_ = sigma;
}

void ElectrostaticFieldSolver::UpdateBoundaryConditions(double t)
{
    const auto &bcs = config_.electric_field.dirichlet_bcs;
    for (size_t i = 0; i < bcs.size(); ++i)
    {
        if (!bcs[i].value_expr.empty())
        {
            double val = bc_expressions_[i].Evaluate({0, 0, 0, t, 0});
            cached_bc_.dirichlet_coeffs[i].constant = val;
        }
    }
}

void ElectrostaticFieldSolver::BuildBCMarkers()
{
    const int num_bdr = config_.fe.GetMesh().bdr_attributes.Max();
    const auto &field_config = config_.electric_field;

    // Dirichlet BCs
    cached_bc_.essential_bdr.SetSize(num_bdr);
    cached_bc_.essential_bdr = 0;
    cached_bc_.dirichlet_coeffs.reserve(field_config.dirichlet_bcs.size());
    cached_bc_.dirichlet_markers.reserve(field_config.dirichlet_bcs.size());
    for (const auto &dbc : field_config.dirichlet_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : dbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
            {
                marker[attr - 1] = 1;
                cached_bc_.essential_bdr[attr - 1] = 1;
            }
        }
        cached_bc_.dirichlet_coeffs.emplace_back(dbc.value);
        cached_bc_.dirichlet_markers.push_back(std::move(marker));
    }

    // Robin BCs
    cached_bc_.robin_l_coeffs.reserve(field_config.robin_bcs.size());
    cached_bc_.robin_q_coeffs.reserve(field_config.robin_bcs.size());
    cached_bc_.robin_markers.reserve(field_config.robin_bcs.size());
    for (const auto &rbc : field_config.robin_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : rbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        cached_bc_.robin_l_coeffs.emplace_back(rbc.l);
        cached_bc_.robin_q_coeffs.emplace_back(rbc.q);
        cached_bc_.robin_markers.push_back(std::move(marker));
    }

    // Cache essential tdofs (won't change)
    config_.fe.GetScalarFESpace().GetEssentialTrueDofs(cached_bc_.essential_bdr,
                                                       cached_essential_tdofs_);
}

ElectrostaticFieldSolver::Coefficients ElectrostaticFieldSolver::BuildCoefficients()
{
    Coefficients c;
    if (sigma_)
    {
        c.sigma = sigma_;
    }
    else
    {
        auto local_sigma = std::make_unique<coeff::ExpressionCoefficient>(
            "electrical_conductivity", config_.materials, config_.fe.GetMesh(), temperature_gf_);
        if (!temperature_gf_)
            local_sigma->SetReferenceTemperature(config_.electric_field.reference_temperature);
        c.sigma = local_sigma.get();
        c.local_sigma = std::move(local_sigma);
    }
    return c;
}

void ElectrostaticFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->debug("Solving electrostatic field");

    auto coeffs = BuildCoefficients();

    auto system = assembly::ElectrostaticAssembler::Assemble(
        config_.fe.GetScalarFESpace(), *coeffs.sigma, coeffs.source, cached_bc_, voltage_);

    if (!cached_solver_)
    {
        cached_solver_ = solver::CreateLinearSolver(config_.simulation.GetSolver("electrostatic"),
                                                    1e-12, 1e-12, 2000, 0);
    }
    cached_solver_->Solve(*system.A, system.B, system.X);
    system.bilinear->RecoverFEMSolution(system.X, *system.linear, voltage_);

    logger->debug("Electrostatic: V=[{:.6e}, {:.6e}]", voltage_.Min(), voltage_.Max());
}

}  // namespace fem::physics

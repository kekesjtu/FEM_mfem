#include "fem/physics/ElectrostaticFieldSolver.hpp"

#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/coeff/Coeffmanagaer.hpp"
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
        {
            bc_expressions_.emplace_back(dbc.value_expr);
        }
        else
        {
            bc_expressions_.emplace_back(std::to_string(dbc.value));
        }
    }
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
    auto &bcs = config_.electric_field.dirichlet_bcs;
    for (size_t i = 0; i < bcs.size(); ++i)
    {
        if (!bcs[i].value_expr.empty())
        {
            bcs[i].value = bc_expressions_[i].Evaluate({0, 0, 0, t, 0});
        }
    }
}

void ElectrostaticFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->info("=== Solving electrostatic field ===");

    const auto &field_config = config_.electric_field;

    // Build electrical conductivity coefficient (sigma), optionally T-dependent.
    // If an external sigma was provided (shared with ThermalFieldSolver for Joule heating),
    // reuse it to avoid evaluating the same expression twice per quadrature point.
    std::unique_ptr<coeff::ExpressionCoefficient> local_sigma;
    if (!sigma_)
    {
        local_sigma = std::make_unique<coeff::ExpressionCoefficient>(
            "electrical_conductivity", config_.materials, config_.fe.GetMesh(), temperature_gf_);
        if (!temperature_gf_)
            local_sigma->SetReferenceTemperature(config_.electric_field.reference_temperature);
        sigma_ = local_sigma.get();
    }

    // Source term (usually 0 for electrostatics: -div(sigma grad V) = 0)
    mfem::ConstantCoefficient source_coeff(0.0);

    // Assemble — boundary marker construction handled by assembler
    assembly::PoissonAssemblyInput input{
        config_.fe.GetScalarFESpace(), *sigma_, source_coeff, field_config.dirichlet_bcs,
        field_config.robin_bcs,        0.0,
    };

    auto system = assembly::MfemPoissonAssembler::Assemble(input, voltage_);

    // Solve
    auto solver = solver::CreateLinearSolver(config_.simulation.GetSolver("electrostatic"), 1e-12,
                                             1e-12, 2000, 0);
    solver->Solve(*system.A, system.B, system.X);

    // Recover the solution
    system.bilinear->RecoverFEMSolution(system.X, *system.linear, voltage_);

    if (!config_.fe.IsParallel())
        logger->info("Electrostatic solve complete. V range: [{:.6e}, {:.6e}]", voltage_.Min(),
                     voltage_.Max());
}

}  // namespace fem::physics

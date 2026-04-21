#include "fem/physics/ElectrostaticFieldSolver.hpp"

#include "fem/assembler/ElectrostaticAssembler.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

#include <stdexcept>

namespace fem::physics
{

ElectrostaticFieldSolver::ElectrostaticFieldSolver(frontend::ProjectConfig &config)
    : config_(config), voltage_(&config.fe.GetScalarFESpace())
{
    voltage_ = 0.0;
    coeffs_.BuildCoeffs(config_, &temperature_gf_, &sigma_);

    // Build and cache BC markers once (they never change)
    cached_bc_ = BCbuilder::ScalarBCBuilder::BuildBC(config_.electric_field, config_.fe.GetMesh(),
                                                     config_.fe.GetScalarFESpace());
}

void ElectrostaticFieldSolver::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    if (!temperature_gf)
        throw std::invalid_argument(
            "ElectrostaticFieldSolver::SetTemperatureField does not accept null");

    temperature_gf_ = temperature_gf;
}

void ElectrostaticFieldSolver::SetElectricalConductivity(mfem::Coefficient *sigma)
{
    if (!sigma)
        throw std::invalid_argument(
            "ElectrostaticFieldSolver::SetElectricalConductivity does not accept null");

    sigma_ = sigma;
}

void ElectrostaticFieldSolver::UpdateBoundaryConditions(double t)
{
    coeffs_.UpdateBoundaryConditions(config_.electric_field, t, cached_bc_);
}

mfem::Coefficient *ElectrostaticFieldSolver::GetSigma()
{
    return coeffs_.GetSigma();
}

const mfem::Coefficient *ElectrostaticFieldSolver::GetSigma() const
{
    return coeffs_.GetSigma();
}

void ElectrostaticFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->debug("Solving electrostatic field");

    auto system = assembler::ElectrostaticAssembler::Assemble(config_.fe.GetScalarFESpace(),
                                                              coeffs_, cached_bc_, voltage_);

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

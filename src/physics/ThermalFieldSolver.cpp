#include "fem/physics/ThermalFieldSolver.hpp"

#include "fem/assembly/MfemPoissonAssembler.hpp"
#include "fem/coeff/Coeffmanagaer.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

namespace fem::physics
{

ThermalFieldSolver::ThermalFieldSolver(frontend::ProjectConfig &config)
    : config_(config), temperature_(&config.fe.GetScalarFESpace())
{
    temperature_ = 293.15;  // initial room temperature
}

void ThermalFieldSolver::SetVoltageField(mfem::GridFunction *voltage_gf)
{
    voltage_gf_ = voltage_gf;
}

void ThermalFieldSolver::SetElectricalConductivity(mfem::Coefficient *sigma)
{
    sigma_ = sigma;
}

void ThermalFieldSolver::EnableTransient(double dt, mfem::GridFunction *T_old)
{
    transient_dt_ = dt;
    T_old_ = T_old;
}

void ThermalFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->info("=== Solving thermal field ===");

    auto &mesh = config_.fe.GetMesh();
    const auto &field_config = config_.thermal_field;

    // Build thermal diffusion coefficient (k)
    coeff::PiecewiseConstantCoefficient k_coeff("diffusion", config_.materials, mesh);

    // Build source coefficient
    std::unique_ptr<mfem::Coefficient> source_ptr;

    if (voltage_gf_ && sigma_)
    {
        // Joule heating: Q = sigma * |grad(V)|^2
        source_ptr = std::make_unique<coeff::JouleHeatingCoefficient>(*sigma_, *voltage_gf_);
        logger->debug("Joule heating source enabled");
    }
    else
    {
        if (!field_config.domain_to_source.empty())
        {
            int max_attr = coeff::GetMaxAttribute(mesh);
            auto pw = std::make_unique<mfem::PWConstCoefficient>(max_attr);
            double default_source = std::stod(field_config.source_default);
            for (int a = 1; a <= max_attr; ++a)
            {
                auto it = field_config.domain_to_source.find(a);
                pw->operator()(a) = (it != field_config.domain_to_source.end())
                                        ? std::stod(it->second)
                                        : default_source;
            }
            source_ptr = std::move(pw);
        }
        else
        {
            source_ptr =
                std::make_unique<mfem::ConstantCoefficient>(std::stod(field_config.source_default));
        }
    }

    // --- Transient backward-Euler mass terms ---
    // LHS adds: (rho_cp/dt * u, v),  RHS adds: (rho_cp/dt * T_old, v)
    std::unique_ptr<coeff::PiecewiseConstantCoefficient> rho_cp_ptr;
    std::unique_ptr<mfem::ConstantCoefficient> inv_dt_ptr;
    std::unique_ptr<mfem::ProductCoefficient> mass_lhs_ptr;
    std::unique_ptr<mfem::GridFunctionCoefficient> T_old_gf_ptr;
    std::unique_ptr<mfem::ProductCoefficient> mass_rhs_ptr;

    mfem::Coefficient *mass_lhs = nullptr;
    mfem::Coefficient *mass_rhs = nullptr;

    if (T_old_ && transient_dt_ > 0.0)
    {
        rho_cp_ptr = std::make_unique<coeff::PiecewiseConstantCoefficient>("rho_cp",
                                                                           config_.materials, mesh);
        inv_dt_ptr = std::make_unique<mfem::ConstantCoefficient>(1.0 / transient_dt_);
        mass_lhs_ptr = std::make_unique<mfem::ProductCoefficient>(*inv_dt_ptr, *rho_cp_ptr);
        T_old_gf_ptr = std::make_unique<mfem::GridFunctionCoefficient>(T_old_);
        mass_rhs_ptr = std::make_unique<mfem::ProductCoefficient>(*mass_lhs_ptr, *T_old_gf_ptr);
        mass_lhs = mass_lhs_ptr.get();
        mass_rhs = mass_rhs_ptr.get();
    }

    // Assemble — boundary marker construction handled by assembler
    assembly::PoissonAssemblyInput input{
        config_.fe.GetScalarFESpace(), k_coeff, *source_ptr, field_config.dirichlet_bcs,
        field_config.robin_bcs,        293.15,  mass_lhs,    mass_rhs,
    };

    auto system = assembly::MfemPoissonAssembler::Assemble(input, temperature_);

    // Solve
    auto solver = solver::CreateLinearSolver(config_.simulation.solver, 1e-12, 1e-12, 2000, 0);
    solver->Solve(*system.A, system.B, system.X);

    // Recover the solution
    system.bilinear->RecoverFEMSolution(system.X, *system.linear, temperature_);

    if (!config_.fe.IsParallel())
        logger->info("Thermal solve complete. T range: [{:.6e}, {:.6e}]", temperature_.Min(),
                     temperature_.Max());
}

}  // namespace fem::physics

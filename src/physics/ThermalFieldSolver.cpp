#include "fem/physics/ThermalFieldSolver.hpp"

#include "fem/assembler/ThermalAssembler.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

#include <stdexcept>

namespace fem::physics
{

ThermalFieldSolver::ThermalFieldSolver(frontend::ProjectConfig &config)
    : config_(config),
      T_np1_(&config.fe.GetScalarFESpace()),
      T_np1_bdf1_(&config.fe.GetScalarFESpace()),
      T_np1_bdf2_(&config.fe.GetScalarFESpace())
{
    T_np1_ = config_.thermal_field.initial_temperature;
    T_np1_bdf1_ = config_.thermal_field.initial_temperature;
    T_np1_bdf2_ = config_.thermal_field.initial_temperature;
    coeffs_.BuildCoeffs(config_, &voltage_, &sigma_);

    cached_bc_ = BCbuilder::ScalarBCBuilder::BuildBC(config_.thermal_field, config_.fe.GetMesh(),
                                                     config_.fe.GetScalarFESpace());
}

void ThermalFieldSolver::SetVoltageField(mfem::GridFunction *voltage_gf)
{
    if (!voltage_gf)
        throw std::invalid_argument("ThermalFieldSolver::SetVoltageField does not accept null");

    voltage_ = voltage_gf;
}

void ThermalFieldSolver::SetElectricalConductivity(mfem::Coefficient *sigma)
{
    if (!sigma)
        throw std::invalid_argument(
            "ThermalFieldSolver::SetElectricalConductivity does not accept null");

    sigma_ = sigma;
}

void ThermalFieldSolver::EnableTransient(double dt, mfem::GridFunction *T_n)
{
    transient_dt_ = dt;
    T_n_ = T_n;
}

void ThermalFieldSolver::ApplyDirichletValues()
{
    T_np1_ = config_.thermal_field.initial_temperature;
    for (size_t i = 0; i < cached_bc_.dirichlet_markers.size(); ++i)
        T_np1_.ProjectBdrCoefficient(cached_bc_.dirichlet_coeffs[i],
                                     cached_bc_.dirichlet_markers[i]);
}

void ThermalFieldSolver::ProjectDirichletBCs()
{
    for (size_t i = 0; i < cached_bc_.dirichlet_markers.size(); ++i)
        T_np1_.ProjectBdrCoefficient(cached_bc_.dirichlet_coeffs[i],
                                     cached_bc_.dirichlet_markers[i]);
}

void ThermalFieldSolver::ApplyBCToRHS(mfem::Vector &b, const mfem::Vector &bc_vec, double alpha_k,
                                      double alpha_c)
{
    K_->Mult(-alpha_k, bc_vec, 1.0, b);
    C_->Mult(-alpha_c, bc_vec, 1.0, b);
    for (int i = 0; i < cached_bc_.essential_tdofs.Size(); i++)
    {
        int td = cached_bc_.essential_tdofs[i];
        b(td) = bc_vec(td);
    }
}

void ThermalFieldSolver::CacheMatrices()
{
    if (matrices_cached_)
        return;

    auto logger = fem::log::Get();
    logger->debug("Caching thermal K and C matrices...");

    auto m = assembler::ThermalAssembler::AssembleMatrices(config_.fe.GetScalarFESpace(), coeffs_,
                                                           cached_bc_);

    K_bf_ = std::move(m.K_bf);
    C_bf_ = std::move(m.C_bf);
    K_ = std::move(m.K);
    C_ = std::move(m.C);

    matrices_cached_ = true;
    logger->debug("Thermal K and C matrices cached.");
}

double ThermalFieldSolver::ComputeInitialRate()
{
    if (!matrices_cached_)
        throw std::runtime_error("CacheMatrices() must be called before ComputeInitialRate()");

    mfem::Vector F_current = assembler::ThermalAssembler::AssembleSourceRHS(
        config_.fe.GetScalarFESpace(), coeffs_, cached_bc_);

    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector T_true(N), residual(N);
    T_np1_.GetTrueDofs(T_true);

    // residual = F - K*T_0
    K_->Mult(-1.0, T_true, 0.0, residual);
    residual += F_current;

    // Approximate dT/dt ≈ diag(C)^{-1} * residual
    mfem::Vector C_diag(N);
    C_->GetDiag(C_diag);

    double local_max = 0.0;
    for (int i = 0; i < N; i++)
    {
        double rate = std::abs(residual(i)) / std::max(C_diag(i), 1e-30);
        local_max = std::max(local_max, rate);
    }

    double global_max = 0.0;
    MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_max;
}

void ThermalFieldSolver::SolveSteady()
{
    auto logger = fem::log::Get();
    logger->debug("Solving thermal field (steady)");

    if (!matrices_cached_)
    {
        logger->critical("CacheMatrices() must be called first");
        throw std::runtime_error("CacheMatrices() must be called first");
    }

    mfem::Vector F_current = assembler::ThermalAssembler::AssembleSourceRHS(
        config_.fe.GetScalarFESpace(), coeffs_, cached_bc_);
    ApplyDirichletValues();

    EnsureSystemMatrix(1.0, 0.0, steady_state_);

    mfem::Vector b(F_current);

    mfem::Vector T_true(config_.fe.GetScalarFESpace().GetTrueVSize());
    ApplyBCAndSolve(*steady_state_.A, b, 1.0, 0.0, steady_state_, &T_true);
    T_np1_.Distribute(T_true);

    logger->debug("Thermal steady: T=[{:.6e}, {:.6e}]", T_np1_.Min(), T_np1_.Max());
}

void ThermalFieldSolver::SolveBDF1()
{
    auto logger = fem::log::Get();
    logger->debug("Solving thermal field (BE)");

    if (!matrices_cached_)
    {
        logger->critical("CacheMatrices() must be called first");
        throw std::runtime_error("CacheMatrices() must be called first");
    }
    if (!T_n_ || transient_dt_ <= 0.0)
    {
        logger->critical("EnableTransient() must be called first");
        throw std::runtime_error("EnableTransient() must be called first");
    }

    const double inv_dt = 1.0 / transient_dt_;

    mfem::Vector F_current = assembler::ThermalAssembler::AssembleSourceRHS(
        config_.fe.GetScalarFESpace(), coeffs_, cached_bc_);
    ApplyDirichletValues();

    EnsureSystemMatrix(1.0, inv_dt, bdf1_state_);

    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector b(N);
    b = 0.0;

    mfem::Vector T_n_true(N);
    T_n_->GetTrueDofs(T_n_true);
    C_->Mult(inv_dt, T_n_true, 0.0, b);
    b += F_current;

    mfem::Vector T_np1_true(N);
    ApplyBCAndSolve(*bdf1_state_.A, b, 1.0, inv_dt, bdf1_state_, &T_np1_true);
    T_np1_.Distribute(T_np1_true);
    T_np1_bdf1_.Distribute(T_np1_true);

    logger->debug("Thermal BE: T=[{:.6e}, {:.6e}]", T_np1_.Min(), T_np1_.Max());
}

void ThermalFieldSolver::SolveBDF2(const mfem::Vector &T_nm1_true, double dt_prev)
{
    auto logger = fem::log::Get();
    logger->debug("Solving thermal field (BDF2)");

    if (!matrices_cached_)
    {
        logger->critical("CacheMatrices() must be called first");
        throw std::runtime_error("CacheMatrices() must be called first");
    }
    if (!T_n_ || transient_dt_ <= 0.0)
    {
        logger->critical("EnableTransient() must be called first");
        throw std::runtime_error("EnableTransient() must be called first");
    }
    if (dt_prev <= 0.0)
    {
        logger->critical("dt_prev must be > 0 for BDF2");
        throw std::runtime_error("dt_prev must be > 0 for BDF2");
    }

    mfem::Vector F_current = assembler::ThermalAssembler::AssembleSourceRHS(
        config_.fe.GetScalarFESpace(), coeffs_, cached_bc_);

    // Only project BCs onto essential DOFs; preserve the Picard iterate in T_np1_ as
    // the initial guess for the linear solver.
    ProjectDirichletBCs();

    // Variable-step BDF2 coefficients
    // r = h_{n+1} / h_n, h = h_{n+1} (current step)
    const double h = transient_dt_;
    const double r = h / dt_prev;
    const double alpha = (1.0 + 2.0 * r) / ((1.0 + r) * h);  // coeff of T_{n+1} in C*dT/dt
    const double beta_n = (1.0 + r) / h;                     // coeff of T_n
    const double beta_nm1 = (r * r) / ((1.0 + r) * h);       // coeff of T_{n-1}

    // System: (K + alpha*C) T_{n+1} = F + C*(beta_n*T_n - beta_nm1*T_{n-1})
    EnsureSystemMatrix(1.0, alpha, bdf2_state_);

    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector b(N);
    b = 0.0;

    mfem::Vector T_n_true(N);
    T_n_->GetTrueDofs(T_n_true);

    C_->Mult(beta_n, T_n_true, 0.0, b);
    C_->Mult(-beta_nm1, T_nm1_true, 1.0, b);
    b += F_current;

    mfem::Vector T_np1_true(N);
    ApplyBCAndSolve(*bdf2_state_.A, b, 1.0, alpha, bdf2_state_, &T_np1_true);
    T_np1_bdf2_.Distribute(T_np1_true);

    logger->debug("Thermal BDF2 solved");
}

void ThermalFieldSolver::EnsureSystemMatrix(double alpha_k, double alpha_c, SolveState &state)
{
    if (state.last_alpha_c != alpha_c || !state.A)
    {
        if (alpha_c == 0.0)
            state.A.reset(new mfem::HypreParMatrix(*K_));
        else
            state.A.reset(mfem::Add(alpha_k, *K_, alpha_c, *C_));
        std::unique_ptr<mfem::HypreParMatrix> elim{
            state.A->EliminateRowsCols(cached_bc_.essential_tdofs)};
        state.solver.reset();
        state.last_alpha_c = alpha_c;
    }
}

void ThermalFieldSolver::ApplyBCAndSolve(mfem::HypreParMatrix &A, mfem::Vector &b, double alpha_k,
                                         double alpha_c, SolveState &state,
                                         mfem::Vector *output_true)
{
    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();

    mfem::Vector x_true(N);
    T_np1_.GetTrueDofs(x_true);

    mfem::Vector bc_vec(N);
    bc_vec = 0.0;
    for (int i = 0; i < cached_bc_.essential_tdofs.Size(); i++)
        bc_vec(cached_bc_.essential_tdofs[i]) = x_true(cached_bc_.essential_tdofs[i]);
    ApplyBCToRHS(b, bc_vec, alpha_k, alpha_c);

    if (!state.solver)
    {
        state.solver = solver::CreateLinearSolver(config_.simulation.GetSolver("thermal"), 1e-12,
                                                  1e-12, 2000, 0);
        state.solver->Solve(A, b, x_true);
    }
    else
    {
        state.solver->Mult(b, x_true);
    }

    *output_true = x_true;
}

}  // namespace fem::physics

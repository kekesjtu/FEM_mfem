#include "fem/physics/ThermalFieldSolver.hpp"

#include "fem/assembly/ThermalAssembler.hpp"
#include "fem/coeff/CoefficientManager.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

namespace fem::physics
{

ThermalFieldSolver::ThermalFieldSolver(frontend::ProjectConfig &config)
    : config_(config), temperature_(&config.fe.GetScalarFESpace())
{
    temperature_ = config_.thermal_field.initial_temperature;

    auto &mesh = config_.fe.GetMesh();
    cached_k_ =
        std::make_unique<coeff::PiecewiseConstantCoefficient>("diffusion", config_.materials, mesh);
    cached_rho_cp_ =
        std::make_unique<coeff::PiecewiseConstantCoefficient>("rho_cp", config_.materials, mesh);

    BuildBCMarkers();
}

void ThermalFieldSolver::SetVoltageField(mfem::GridFunction *voltage_gf)
{
    voltage_ = voltage_gf;
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

void ThermalFieldSolver::BuildBCMarkers()
{
    const int num_bdr = config_.fe.GetMesh().bdr_attributes.Max();
    const auto &field_config = config_.thermal_field;

    cached_essential_bdr_.SetSize(num_bdr);
    cached_essential_bdr_ = 0;
    cached_dirichlet_coeffs_.reserve(field_config.dirichlet_bcs.size());
    cached_dirichlet_markers_.reserve(field_config.dirichlet_bcs.size());
    for (const auto &dbc : field_config.dirichlet_bcs)
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
        cached_dirichlet_coeffs_.emplace_back(dbc.value);
        cached_dirichlet_markers_.push_back(std::move(marker));
    }

    cached_robin_l_coeffs_.reserve(field_config.robin_bcs.size());
    cached_robin_q_coeffs_.reserve(field_config.robin_bcs.size());
    cached_robin_markers_.reserve(field_config.robin_bcs.size());
    for (const auto &rbc : field_config.robin_bcs)
    {
        mfem::Array<int> marker(num_bdr);
        marker = 0;
        for (int attr : rbc.bdr_attributes)
        {
            if (attr >= 1 && attr <= num_bdr)
                marker[attr - 1] = 1;
        }
        cached_robin_l_coeffs_.emplace_back(rbc.l);
        cached_robin_q_coeffs_.emplace_back(rbc.q);
        cached_robin_markers_.push_back(std::move(marker));
    }

    config_.fe.GetScalarFESpace().GetEssentialTrueDofs(cached_essential_bdr_,
                                                       cached_essential_tdofs_);
}

std::unique_ptr<mfem::Coefficient> ThermalFieldSolver::BuildSourceCoefficient()
{
    if (voltage_ && sigma_)
    {
        return std::make_unique<coeff::JouleHeatingCoefficient>(*sigma_, *voltage_);
    }

    const auto &field_config = config_.thermal_field;
    if (!field_config.domain_to_source.empty())
    {
        int max_attr = coeff::GetMaxAttribute(config_.fe.GetMesh());
        auto pw = std::make_unique<mfem::PWConstCoefficient>(max_attr);
        double default_source = std::stod(field_config.source_default);
        for (int a = 1; a <= max_attr; ++a)
        {
            auto it = field_config.domain_to_source.find(a);
            pw->operator()(a) = (it != field_config.domain_to_source.end()) ? std::stod(it->second)
                                                                            : default_source;
        }
        return pw;
    }

    return std::make_unique<mfem::ConstantCoefficient>(std::stod(field_config.source_default));
}

void ThermalFieldSolver::ApplyDirichletValues()
{
    temperature_ = config_.thermal_field.initial_temperature;
    for (size_t i = 0; i < cached_dirichlet_markers_.size(); ++i)
        temperature_.ProjectBdrCoefficient(cached_dirichlet_coeffs_[i],
                                           cached_dirichlet_markers_[i]);
}

void ThermalFieldSolver::ApplyBCToRHS(mfem::Vector &b, const mfem::Vector &bc_vec, double alpha_k,
                                      double alpha_c)
{
    K_->Mult(-alpha_k, bc_vec, 1.0, b);
    C_->Mult(-alpha_c, bc_vec, 1.0, b);
    for (int i = 0; i < cached_essential_tdofs_.Size(); i++)
    {
        int td = cached_essential_tdofs_[i];
        b(td) = bc_vec(td);
    }
}

void ThermalFieldSolver::CacheMatrices()
{
    if (matrices_cached_)
        return;

    auto logger = fem::log::Get();
    logger->debug("Caching thermal K and C matrices...");

    auto m = assembly::ThermalAssembler::AssembleMatrices(config_.fe.GetScalarFESpace(), *cached_k_,
                                                          *cached_rho_cp_, cached_robin_l_coeffs_,
                                                          cached_robin_markers_);

    K_bf_ = std::move(m.K_bf);
    C_bf_ = std::move(m.C_bf);
    K_ = std::move(m.K);
    C_ = std::move(m.C);

    matrices_cached_ = true;
    logger->debug("Thermal K and C matrices cached.");
}

void ThermalFieldSolver::AssembleSourceRHS()
{
    auto source = BuildSourceCoefficient();

    F_current_ = assembly::ThermalAssembler::AssembleSourceRHS(
        config_.fe.GetScalarFESpace(), *source, cached_robin_q_coeffs_, cached_robin_markers_);
}

double ThermalFieldSolver::ComputeInitialRate()
{
    if (!matrices_cached_)
        throw std::runtime_error("CacheMatrices() must be called before ComputeInitialRate()");

    AssembleSourceRHS();

    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector T_true(N), residual(N);
    temperature_.GetTrueDofs(T_true);

    // residual = F - K*T_0
    K_->Mult(-1.0, T_true, 0.0, residual);
    residual += F_current_;

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

    AssembleSourceRHS();
    ApplyDirichletValues();

    EnsureSystemMatrix(1.0, 0.0, steady_state_);

    mfem::Vector b(F_current_);

    mfem::Vector T_true(config_.fe.GetScalarFESpace().GetTrueVSize());
    ApplyBCAndSolve(*steady_state_.A, b, 1.0, 0.0, steady_state_, &T_true);
    temperature_.Distribute(T_true);

    logger->debug("Thermal steady: T=[{:.6e}, {:.6e}]", temperature_.Min(), temperature_.Max());
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
    if (!T_old_ || transient_dt_ <= 0.0)
    {
        logger->critical("EnableTransient() must be called first");
        throw std::runtime_error("EnableTransient() must be called first");
    }

    const double inv_dt = 1.0 / transient_dt_;

    AssembleSourceRHS();
    ApplyDirichletValues();

    EnsureSystemMatrix(1.0, inv_dt, bdf1_state_);

    const int N = config_.fe.GetScalarFESpace().GetTrueVSize();
    mfem::Vector b(N);
    b = 0.0;

    mfem::Vector T_old_true(N);
    T_old_->GetTrueDofs(T_old_true);
    C_->Mult(inv_dt, T_old_true, 0.0, b);
    b += F_current_;

    mfem::Vector T_true(N);
    ApplyBCAndSolve(*bdf1_state_.A, b, 1.0, inv_dt, bdf1_state_, &T_true);
    temperature_.Distribute(T_true);

    logger->debug("Thermal BE: T=[{:.6e}, {:.6e}]", temperature_.Min(), temperature_.Max());
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
    if (!T_old_ || transient_dt_ <= 0.0)
    {
        logger->critical("EnableTransient() must be called first");
        throw std::runtime_error("EnableTransient() must be called first");
    }
    if (dt_prev <= 0.0)
    {
        logger->critical("dt_prev must be > 0 for BDF2");
        throw std::runtime_error("dt_prev must be > 0 for BDF2");
    }

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

    mfem::Vector T_old_true(N);
    T_old_->GetTrueDofs(T_old_true);

    C_->Mult(beta_n, T_old_true, 0.0, b);
    C_->Mult(-beta_nm1, T_nm1_true, 1.0, b);
    b += F_current_;

    T_bdf2_true_.SetSize(N);
    ApplyBCAndSolve(*bdf2_state_.A, b, 1.0, alpha, bdf2_state_, &T_bdf2_true_);

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
        delete state.A->EliminateRowsCols(cached_essential_tdofs_);
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
    temperature_.GetTrueDofs(x_true);

    mfem::Vector bc_vec(N);
    bc_vec = 0.0;
    for (int i = 0; i < cached_essential_tdofs_.Size(); i++)
        bc_vec(cached_essential_tdofs_[i]) = x_true(cached_essential_tdofs_[i]);
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

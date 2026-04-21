#pragma once

#include "fem/BCbuilder/ScalarBCBuilder.hpp"
#include "fem/coeff/ThermalCoeffs.hpp"
#include "fem/frontend/Config.hpp"
#include "fem/solver/ILinearSolver.hpp"
#include "mfem.hpp"

#include <memory>

namespace fem::physics
{

class ThermalFieldSolver
{
  public:
    explicit ThermalFieldSolver(frontend::ProjectConfig &config);

    void SetVoltageField(mfem::GridFunction *voltage_gf);
    void SetElectricalConductivity(mfem::Coefficient *sigma);
    void EnableTransient(double dt, mfem::GridFunction *T_n);
    void CacheMatrices();
    void SolveSteady();
    void SolveBDF1();
    void SolveBDF2(const mfem::Vector &T_nm1_true, double dt_prev);

    mfem::ParGridFunction &GetTemperature()
    {
        return T_np1_;
    }
    mfem::ParGridFunction &GetBDF1Solution()
    {
        return T_np1_bdf1_;
    }
    mfem::ParGridFunction &GetBDF2Solution()
    {
        return T_np1_bdf2_;
    }
    double ComputeInitialRate();

  private:
    void ApplyDirichletValues();
    void ProjectDirichletBCs();
    void ApplyBCToRHS(mfem::Vector &b, const mfem::Vector &bc_vec, double alpha_k, double alpha_c);

    struct SolveState
    {
        std::unique_ptr<mfem::HypreParMatrix> A;
        std::unique_ptr<solver::ILinearSolver> solver;
        double last_alpha_c = -1.0;
    };

    void EnsureSystemMatrix(double alpha_k, double alpha_c, SolveState &state);
    void ApplyBCAndSolve(mfem::HypreParMatrix &A, mfem::Vector &b, double alpha_k, double alpha_c,
                         SolveState &state, mfem::Vector *output_true);

    frontend::ProjectConfig &config_;

    mfem::GridFunction *voltage_ = nullptr;
    mfem::Coefficient *sigma_ = nullptr;
    double transient_dt_ = 0.0;
    mfem::GridFunction *T_n_ = nullptr;
    coeff::ThermalCoeffs coeffs_;

    BCbuilder::ScalarBCData cached_bc_;

    bool matrices_cached_ = false;
    std::unique_ptr<mfem::ParBilinearForm> K_bf_;
    std::unique_ptr<mfem::ParBilinearForm> C_bf_;
    std::unique_ptr<mfem::HypreParMatrix> K_;
    std::unique_ptr<mfem::HypreParMatrix> C_;

    mfem::ParGridFunction T_np1_;
    mfem::ParGridFunction T_np1_bdf1_;
    mfem::ParGridFunction T_np1_bdf2_;

    SolveState steady_state_;
    SolveState bdf1_state_;
    SolveState bdf2_state_;
};

}  // namespace fem::physics

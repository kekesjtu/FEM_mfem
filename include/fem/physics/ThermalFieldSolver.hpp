#pragma once

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
    void EnableTransient(double dt, mfem::GridFunction *T_old);
    void CacheMatrices();
    void SolveSteady();
    void SolveBDF1();
    void SolveBDF2(const mfem::Vector &T_nm1_true, double dt_prev);
    void AssembleSourceRHS();

    mfem::ParGridFunction &GetTemperature()
    {
        return temperature_;
    }
    const mfem::Vector &GetBDF2Solution() const
    {
        return T_bdf2_true_;
    }
    double ComputeInitialRate();

  private:
    void BuildBCMarkers();
    std::unique_ptr<mfem::Coefficient> BuildSourceCoefficient();
    void ApplyDirichletValues();
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
    mfem::ParGridFunction temperature_;
    mfem::GridFunction *voltage_ = nullptr;
    mfem::Coefficient *sigma_ = nullptr;
    double transient_dt_ = 0.0;
    mfem::GridFunction *T_old_ = nullptr;

    std::unique_ptr<mfem::Coefficient> cached_k_;
    std::unique_ptr<mfem::Coefficient> cached_rho_cp_;

    mfem::Array<int> cached_essential_bdr_;
    mfem::Array<int> cached_essential_tdofs_;
    std::vector<mfem::ConstantCoefficient> cached_dirichlet_coeffs_;
    std::vector<mfem::Array<int>> cached_dirichlet_markers_;
    std::vector<mfem::ConstantCoefficient> cached_robin_l_coeffs_;
    std::vector<mfem::ConstantCoefficient> cached_robin_q_coeffs_;
    std::vector<mfem::Array<int>> cached_robin_markers_;

    bool matrices_cached_ = false;
    std::unique_ptr<mfem::ParBilinearForm> K_bf_;
    std::unique_ptr<mfem::ParBilinearForm> C_bf_;
    std::unique_ptr<mfem::HypreParMatrix> K_;
    std::unique_ptr<mfem::HypreParMatrix> C_;

    mfem::Vector F_current_;
    mfem::Vector T_bdf2_true_;

    SolveState steady_state_;
    SolveState bdf1_state_;
    SolveState bdf2_state_;
};

}  // namespace fem::physics

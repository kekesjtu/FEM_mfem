#include "fem/physics/MechanicalFieldSolver.hpp"

#include "fem/assembly/MfemLinearElasticityAssembler.hpp"
#include "fem/coeff/Coeffmanagaer.hpp"
#include "fem/log/Logger.hpp"
#include "fem/solver/LinearSolverFactory.hpp"

namespace fem::physics
{

MechanicalFieldSolver::MechanicalFieldSolver(frontend::ProjectConfig &config)
    : config_(config), displacement_(&config.fe.GetVectorFESpace())
{
    displacement_ = 0.0;
}

void MechanicalFieldSolver::SetTemperatureField(mfem::GridFunction *temperature_gf)
{
    temperature_gf_ = temperature_gf;
}

void MechanicalFieldSolver::Solve()
{
    auto logger = fem::log::Get();
    logger->info("=== Solving mechanical field ===");

    auto &mesh = config_.fe.GetMesh();
    const auto &field_config = config_.mechanical_field;
    const int dim = mesh.Dimension();

    // Lame parameters — auto-computed by ConfigLoader from young_modulus + poisson_ratio
    coeff::PiecewiseConstantCoefficient lambda_coeff("lambda_lame", config_.materials, mesh);
    coeff::PiecewiseConstantCoefficient mu_coeff("mu_lame", config_.materials, mesh);

    // Body force — domain-aware (default + per-domain override, like thermal source)
    std::unique_ptr<mfem::VectorCoefficient> body_force_ptr;
    if (!field_config.domain_to_body_force.empty())
    {
        int max_attr = coeff::GetMaxAttribute(mesh);
        // Build per-component PWConstCoefficient arrays
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
                {
                    comp_vals(a) = it->second[d];
                }
                else
                {
                    comp_vals(a) = default_val;
                }
            }
            pw_vec->Set(d, new mfem::PWConstCoefficient(comp_vals));
        }
        body_force_ptr = std::move(pw_vec);
    }
    else
    {
        mfem::Vector body_force_vec(dim);
        for (int d = 0; d < dim && d < static_cast<int>(field_config.body_force_default.size());
             ++d)
        {
            body_force_vec(d) = field_config.body_force_default[d];
        }
        body_force_ptr = std::make_unique<mfem::VectorConstantCoefficient>(body_force_vec);
    }

    // Optional thermal expansion coupling
    std::unique_ptr<coeff::PiecewiseConstantCoefficient> alpha_coeff;
    std::unique_ptr<mfem::GridFunctionCoefficient> temp_coeff;

    // Assemble — boundary marker construction handled by assembler
    assembly::LinearElasticityInput input{
        config_.fe.GetVectorFESpace(),
        lambda_coeff,
        mu_coeff,
        *body_force_ptr,
        nullptr,  // thermal_expansion_secant
        nullptr,  // temperature
        field_config.reference_temperature,
        field_config.displacement_bcs,
        field_config.normal_displacement_bcs,
        field_config.pressure_bcs,
    };

    if (temperature_gf_)
    {
        alpha_coeff = std::make_unique<coeff::PiecewiseConstantCoefficient>(
            "thermal_expansion_secant", config_.materials, mesh);
        temp_coeff = std::make_unique<mfem::GridFunctionCoefficient>(temperature_gf_);

        input.thermal_expansion_secant = alpha_coeff.get();
        input.temperature = temp_coeff.get();

        logger->debug("Thermal-mechanical coupling enabled, T_ref={:.2f}",
                      field_config.reference_temperature);
    }

    auto system = assembly::MfemLinearElasticityAssembler::Assemble(input, displacement_);

    // Solve
    auto solver = solver::CreateLinearSolver(config_.simulation.solver, 1e-12, 1e-12, 5000, 0);
    solver->Solve(*system.A, system.B, system.X);

    // Recover the solution
    system.bilinear->RecoverFEMSolution(system.X, *system.linear, displacement_);

    if (!config_.fe.IsParallel())
        logger->info("Mechanical solve complete. U range: [{:.6e}, {:.6e}]", displacement_.Min(),
                     displacement_.Max());
}

}  // namespace fem::physics

#include "fem/frontend/ConfigLoader.hpp"

#include "fem/frontend/Expression.hpp"
#include "fem/log/Logger.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

namespace fem::frontend
{
using json = nlohmann::json;

namespace
{
template <typename T>
T ReadOrDefault(const json &node, const std::string &key, T fallback)
{
    if (!node.contains(key))
    {
        return fallback;
    }
    return node.at(key).get<T>();
}
}  // namespace

ProjectConfig ConfigLoader::LoadFromFile(const std::string &path)
{
    std::ifstream ifs(path);
    if (!ifs.is_open())
    {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    json root = json::parse(ifs);
    ConfigLoader loader;
    loader.LoadSimulation(root);

    // Initialize logger early so subsequent steps can log
    fem::log::Init(loader.config_.simulation.log_level);
    auto logger = fem::log::Get();
    logger->info("Loading config from: {}", path);

    loader.LoadMaterials(root);
    loader.LoadFields(root);
    loader.BuildFE();

    logger->info("Config loaded. Electrostatic={}, Thermal={}, Mechanical={}",
                 loader.config_.HasElectricField(), loader.config_.HasThermalField(),
                 loader.config_.HasMechanicalField());
    logger->info("Solver for each field: Electrostatic='{}', Thermal='{}', Mechanical='{}'",
                 loader.config_.simulation.solver_electrostatic,
                 loader.config_.simulation.solver_thermal,
                 loader.config_.simulation.solver_mechanical);
    logger->info("Is nonlinear : {}", loader.config_.NeedsPicardIteration() ? "true" : "false");
    logger->info("Is transient problem: {}", loader.config_.simulation.transient_enabled);

    return std::move(loader.config_);
}

void ConfigLoader::LoadSimulation(const json &root)
{
    if (!root.contains("simulation"))
    {
        throw std::runtime_error("Missing 'simulation' section in config.");
    }

    const auto &sim = root.at("simulation");
    auto &s = config_.simulation;
    s.mesh_path = sim.at("mesh_path").get<std::string>();
    s.order = ReadOrDefault(sim, "order", 1);
    s.uniform_refinement_levels = ReadOrDefault(sim, "uniform_refinement_levels", 0);
    s.output_dir = ReadOrDefault<std::string>(sim, "output_dir", "results");
    s.log_level = ReadOrDefault<std::string>(sim, "log_level", "info");
    s.solver = ReadOrDefault<std::string>(sim, "solver", "pcg");
    s.solver_electrostatic = ReadOrDefault<std::string>(sim, "solver_electrostatic", s.solver);
    s.solver_thermal = ReadOrDefault<std::string>(sim, "solver_thermal", s.solver);
    s.solver_mechanical = ReadOrDefault<std::string>(sim, "solver_mechanical", s.solver);
    s.picard_max_iterations = ReadOrDefault(sim, "picard_max_iterations", 30);
    s.picard_tolerance = ReadOrDefault(sim, "picard_tolerance", 1.0e-6);
    s.picard_relaxation = ReadOrDefault(sim, "picard_relaxation", 1.0);
    s.transient_enabled = ReadOrDefault(sim, "transient_enabled", false);
    s.transient_t_start = ReadOrDefault(sim, "transient_t_start", 0.0);
    s.transient_t_end = ReadOrDefault(sim, "transient_t_end", 1.0);
    s.transient_dt = ReadOrDefault(sim, "transient_dt", 1.0);
    s.transient_output_interval = ReadOrDefault(sim, "transient_output_interval", 10.0);
    s.adaptive_dt = ReadOrDefault(sim, "adaptive_dt", false);
    s.adaptive_reltol = ReadOrDefault(sim, "adaptive_reltol", 1.0e-3);
    s.adaptive_abstol = ReadOrDefault(sim, "adaptive_abstol", 1.0e-6);
    s.dt_min = ReadOrDefault(sim, "dt_min", 1.0e-6);
    s.dt_max = ReadOrDefault(sim, "dt_max", 0.0);
    s.eta_safety = ReadOrDefault(sim, "eta_safety", 0.9);
    s.eta_max = ReadOrDefault(sim, "eta_max", 2.0);
    s.eta_min = ReadOrDefault(sim, "eta_min", 0.2);
    s.comsol_reference_path = ReadOrDefault<std::string>(sim, "comsol_reference_path", "");
    s.compare_args = ReadOrDefault<std::string>(sim, "compare_args", "");

    // --- Adaptive parameter validation ---
    if (s.transient_enabled)
    {
        auto logger = fem::log::Get();
        if (s.adaptive_reltol <= 0.0 && s.adaptive_abstol <= 0.0)
            throw std::runtime_error("adaptive_reltol and adaptive_abstol cannot both be <= 0");
        if (s.adaptive_reltol < 0.0)
            throw std::runtime_error("adaptive_reltol must be >= 0");
        if (s.adaptive_abstol < 0.0)
            throw std::runtime_error("adaptive_abstol must be >= 0");
        if (s.dt_min <= 0.0)
            throw std::runtime_error("dt_min must be > 0");
        double eff_dt_max = (s.dt_max > 0.0) ? s.dt_max : s.transient_output_interval;
        if (s.dt_min > eff_dt_max)
            throw std::runtime_error("dt_min (" + std::to_string(s.dt_min) +
                                     ") must be <= dt_max (" + std::to_string(eff_dt_max) + ")");
        if (s.eta_max <= 1.0)
            throw std::runtime_error("eta_max must be > 1.0");
        if (s.eta_min <= 0.0 || s.eta_min >= 1.0)
            throw std::runtime_error("eta_min must be in (0, 1)");
        if (s.eta_safety <= 0.0 || s.eta_safety > 1.0)
            throw std::runtime_error("eta_safety must be in (0, 1]");
        if (s.eta_min >= s.eta_max)
            throw std::runtime_error("eta_min must be < eta_max");
    }
}

void ConfigLoader::LoadMaterials(const json &root)
{
    if (!root.contains("materials") || !root.contains("domain_materials"))
    {
        return;
    }

    auto &db = config_.materials;
    const auto &mats = root.at("materials");
    for (const auto &mat : mats)
    {
        const std::string name = mat.at("name").get<std::string>();
        const auto &props = mat.at("properties");
        for (auto it = props.begin(); it != props.end(); ++it)
        {
            db.material_to_properties[name][it.key()] = it.value().get<std::string>();
        }
    }

    const auto &dm = root.at("domain_materials");
    for (auto it = dm.begin(); it != dm.end(); ++it)
    {
        const std::string &mat_name = it.key();
        std::vector<int> domains = it.value().get<std::vector<int>>();
        db.domain_to_material[mat_name] = std::move(domains);
    }

    // Auto-compute Lame parameters from Young's modulus and Poisson ratio
    for (auto &[mat_name, props] : db.material_to_properties)
    {
        auto E_it = props.find("young_modulus");
        auto nu_it = props.find("poisson_ratio");
        if (E_it != props.end() && nu_it != props.end())
        {
            double E = std::stod(E_it->second);
            double nu = std::stod(nu_it->second);
            double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
            double mu = E / (2.0 * (1.0 + nu));
            props["lambda_lame"] = std::to_string(lambda);
            props["mu_lame"] = std::to_string(mu);
        }

        // Auto-compute rho_cp from density and heat_capacity
        auto rho_it = props.find("density");
        auto cp_it = props.find("heat_capacity");
        if (rho_it != props.end() && cp_it != props.end())
        {
            double rho = std::stod(rho_it->second);
            double cp = std::stod(cp_it->second);
            props["rho_cp"] = std::to_string(rho * cp);
        }
    }
}

void ConfigLoader::LoadFields(const json &root)
{
    if (!root.contains("fields"))
    {
        throw std::runtime_error("Missing 'fields' section in config.");
    }

    for (const auto &field : root.at("fields"))
    {
        const std::string type = field.at("type").get<std::string>();
        if (type == "electrostatic")
        {
            config_.electric_field = LoadScalarField(field);
            config_.electric_field.type = "electrostatic";
            config_.electric_field.enabled = true;
        }
        else if (type == "thermal")
        {
            config_.thermal_field = LoadScalarField(field);
            config_.thermal_field.type = "thermal";
            config_.thermal_field.enabled = true;
        }
        else if (type == "mechanical")
        {
            config_.mechanical_field = LoadMechanicalField(field);
            config_.mechanical_field.type = "mechanical";
            config_.mechanical_field.enabled = true;
        }
        else
        {
            throw std::runtime_error("Unknown field type: " + type);
        }
    }

    // Propagate thermal initial_temperature → electric/mechanical reference_temperature
    if (config_.thermal_field.enabled)
    {
        const double T_init = config_.thermal_field.initial_temperature;
        config_.electric_field.reference_temperature = T_init;
        config_.mechanical_field.reference_temperature = T_init;
    }
}

ScalarFieldConfig ConfigLoader::LoadScalarField(const json &field_node)
{
    ScalarFieldConfig cfg;
    cfg.source_default = ReadOrDefault<std::string>(field_node, "source_default", "0.0");

    if (field_node.contains("source_by_domain"))
    {
        for (auto it = field_node.at("source_by_domain").begin();
             it != field_node.at("source_by_domain").end(); ++it)
        {
            int domain_id = std::stoi(it.key());
            cfg.domain_to_source[domain_id] = it.value().get<std::string>();
        }
    }

    if (field_node.contains("dirichlet_boundary_conditions"))
    {
        for (const auto &bc : field_node.at("dirichlet_boundary_conditions"))
        {
            ScalarFieldConfig::DirichletBC dbc;
            dbc.bdr_attributes = bc.at("bdr_attributes").get<std::vector<int>>();
            const auto &val = bc.at("value");
            if (val.is_string())
            {
                dbc.value_expr = val.get<std::string>();
                // Evaluate at t=0 for initial value
                frontend::Expression expr(dbc.value_expr);
                dbc.value = expr.Evaluate({0, 0, 0, 0, 0});
            }
            else
            {
                dbc.value = val.get<double>();
            }
            cfg.dirichlet_bcs.push_back(std::move(dbc));
        }
    }

    if (field_node.contains("robin_boundary_conditions"))
    {
        for (const auto &bc : field_node.at("robin_boundary_conditions"))
        {
            ScalarFieldConfig::RobinBC rbc;
            rbc.bdr_attributes = bc.at("bdr_attributes").get<std::vector<int>>();
            rbc.l = bc.at("l").get<double>();
            rbc.q = bc.at("q").get<double>();
            cfg.robin_bcs.push_back(std::move(rbc));
        }
    }

    cfg.reference_temperature = ReadOrDefault(field_node, "reference_temperature", 293.15);
    cfg.initial_temperature = ReadOrDefault(field_node, "initial_temperature", 293.15);

    return cfg;
}

MechanicalFieldConfig ConfigLoader::LoadMechanicalField(const json &field_node)
{
    MechanicalFieldConfig cfg;

    if (field_node.contains("body_force_default"))
    {
        cfg.body_force_default = field_node.at("body_force_default").get<std::vector<double>>();
    }

    if (field_node.contains("body_force_by_domain"))
    {
        for (auto it = field_node.at("body_force_by_domain").begin();
             it != field_node.at("body_force_by_domain").end(); ++it)
        {
            int domain_id = std::stoi(it.key());
            cfg.domain_to_body_force[domain_id] = it.value().get<std::vector<double>>();
        }
    }

    if (field_node.contains("displacement_boundary_conditions"))
    {
        for (const auto &bc : field_node.at("displacement_boundary_conditions"))
        {
            MechanicalFieldConfig::DisplacementBC dbc;
            dbc.bdr_attributes = bc.at("bdr_attributes").get<std::vector<int>>();
            dbc.value = bc.at("value").get<std::vector<double>>();
            cfg.displacement_bcs.push_back(std::move(dbc));
        }
    }

    if (field_node.contains("normal_displacement_boundary_conditions"))
    {
        for (const auto &bc : field_node.at("normal_displacement_boundary_conditions"))
        {
            MechanicalFieldConfig::NormalDisplacementBC nbc;
            nbc.bdr_attributes = bc.at("bdr_attributes").get<std::vector<int>>();
            nbc.penalty = ReadOrDefault(bc, "penalty", 1.0e12);
            cfg.normal_displacement_bcs.push_back(std::move(nbc));
        }
    }

    if (field_node.contains("pressure_boundary_conditions"))
    {
        for (const auto &bc : field_node.at("pressure_boundary_conditions"))
        {
            MechanicalFieldConfig::PressureBC pbc;
            pbc.bdr_attributes = bc.at("bdr_attributes").get<std::vector<int>>();
            pbc.value = bc.at("value").get<double>();
            cfg.pressure_bcs.push_back(std::move(pbc));
        }
    }

    cfg.reference_temperature = ReadOrDefault(field_node, "reference_temperature", 293.15);

    return cfg;
}

void ConfigLoader::BuildFE()
{
    auto logger = fem::log::Get();
    const auto &sim = config_.simulation;

    // Load serial mesh first
    logger->info("Loading mesh from: {}", sim.mesh_path);
    config_.fe.serial_mesh = std::make_unique<mfem::Mesh>(sim.mesh_path.c_str(), 1, 1);

    for (int i = 0; i < sim.uniform_refinement_levels; ++i)
    {
        config_.fe.serial_mesh->UniformRefinement();
    }

    const int dim = config_.fe.serial_mesh->Dimension();
    logger->info("Mesh loaded: dim={}, elements={}, vertices={}, bdr_elements={}", dim,
                 config_.fe.serial_mesh->GetNE(), config_.fe.serial_mesh->GetNV(),
                 config_.fe.serial_mesh->GetNBE());

    // Always create ParMesh (works with np=1 too)
    logger->info("Execution mode: parallel (MPI ranks={})", mfem::Mpi::WorldSize());
    config_.fe.pmesh = std::make_unique<mfem::ParMesh>(MPI_COMM_WORLD, *config_.fe.serial_mesh);
    logger->info("ParMesh created: local elements={}", config_.fe.pmesh->GetNE());

    if (config_.HasElectricField() || config_.HasThermalField())
    {
        config_.fe.scalar_fec = std::make_unique<mfem::H1_FECollection>(sim.order, dim);
        config_.fe.scalar_fespace = std::make_unique<mfem::ParFiniteElementSpace>(
            config_.fe.pmesh.get(), config_.fe.scalar_fec.get());
        logger->info("Scalar FE space: order={}, global_ndof={}", sim.order,
                     config_.fe.scalar_fespace->GlobalTrueVSize());
    }

    if (config_.HasMechanicalField())
    {
        config_.fe.vector_fec = std::make_unique<mfem::H1_FECollection>(sim.order, dim);
        config_.fe.vector_fespace = std::make_unique<mfem::ParFiniteElementSpace>(
            config_.fe.pmesh.get(), config_.fe.vector_fec.get(), dim);
        logger->info("Vector FE space: order={}, global_ndof={}", sim.order,
                     config_.fe.vector_fespace->GlobalTrueVSize());
    }
}

bool ProjectConfig::NeedsPicardIteration() const
{
    // Picard iteration is needed when electrostatic and thermal are both enabled AND
    // electrical_conductivity depends on T (i.e. is not a pure constant)
    if (!HasElectricField() || !HasThermalField())
    {
        return false;
    }

    for (const auto &[mat_name, props] : materials.material_to_properties)
    {
        auto it = props.find("electrical_conductivity");
        if (it != props.end())
        {
            const std::string &expr = it->second;
            frontend::Expression test_expr(expr);
            if (test_expr.UsesVariable("T"))
            {
                return true;
            }
        }
    }
    return false;
}
}  // namespace fem::frontend

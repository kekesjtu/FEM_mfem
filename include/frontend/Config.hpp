#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include "mfem.hpp"

namespace fem::frontend
{

// 仿真参数
struct SimulationConfig
{
    std::string mesh_path;
    int order = 1;
    int uniform_refinement_levels = 0;
    std::string output_dir = "results";
    std::string log_level = "info";
    std::string solver = "pcg";
    std::string solver_electrostatic;  // per-field override, empty → use solver
    std::string solver_thermal;        // per-field override, empty → use solver
    std::string solver_mechanical;     // per-field override, empty → use solver
    int picard_max_iterations = 30;
    double picard_tolerance = 1.0e-6;
    double picard_relaxation = 1.0;
    bool transient_enabled = false;
    double transient_t_start = 0.0;
    double transient_t_end = 1.0;
    double transient_dt = 1.0;
    double transient_output_interval = 10.0;
    bool adaptive_dt = false;
    double adaptive_reltol = 1.0e-3;
    double adaptive_abstol = 1.0e-6;
    double dt_min = 1.0e-6;
    double dt_max = 0.0;  // 0 = use transient_output_interval
    double eta_safety = 0.9;
    double eta_max = 2.0;
    double eta_min = 0.2;
    std::string comsol_reference_path;
    std::string compare_args;
    std::string curve_reference_path;
    std::string curve_compare_args;

    /// Return the effective solver name for a given field type.
    const std::string &GetSolver(const std::string &field_type) const
    {
        if (field_type == "electrostatic" && !solver_electrostatic.empty())
            return solver_electrostatic;
        if (field_type == "thermal" && !solver_thermal.empty())
            return solver_thermal;
        if (field_type == "mechanical" && !solver_mechanical.empty())
            return solver_mechanical;
        return solver;
    }
};

// 材料数据结构
struct MaterialDatabase
{
    // 材料名 -> MFEM volume attribute id 列表。
    std::unordered_map<std::string, std::vector<int>> material_to_attributes;

    // 材料名 -> 材料属性表达式表。
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
        material_to_properties;

    // LoadMaterials() 构造的 attribute id 索引，后续系数代码不再重复构造。
    std::unordered_map<int, std::string> attr_to_material;
    std::unordered_map<int, std::unordered_map<std::string, std::string>> attr_to_properties;
};

struct ScalarFieldConfig
{
    bool enabled = false;

    struct DirichletBC
    {
        std::vector<int> bdr_attributes;
        double value = 0.0;
        /// Optional muparser expression for time-varying BCs (e.g. "0.02*sin(100*pi*t)").
        /// If non-empty, value is re-evaluated at each transient time step.
        std::string value_expr;
    };

    struct RobinBC
    {
        std::vector<int> bdr_attributes;
        double l = 0.0;
        double q = 0.0;
    };

    std::string type;
    std::string source_default = "0.0";
    std::unordered_map<int, std::string> attr_to_source;
    std::vector<DirichletBC> dirichlet_bcs;
    std::vector<RobinBC> robin_bcs;
    /// For electrostatic fields: temperature used when evaluating T-dependent conductivity
    /// expressions in single-field mode or as the starting temperature before Picard coupling.
    double reference_temperature = 293.15;
    /// For thermal fields: initial temperature (doubles as reference for electrostatic/mechanical).
    double initial_temperature = 293.15;
};

struct MechanicalFieldConfig
{
    bool enabled = false;

    struct DisplacementBC
    {
        std::vector<int> bdr_attributes;
        std::vector<double> value{0.0, 0.0, 0.0};
    };

    struct NormalDisplacementBC
    {
        std::vector<int> bdr_attributes;
        double penalty = 1.0e12;
    };

    struct PressureBC
    {
        std::vector<int> bdr_attributes;
        double value = 0.0;
    };

    std::string type;
    std::vector<double> body_force_default{0.0, 0.0, 0.0};
    std::unordered_map<int, std::vector<double>> attr_to_body_force;
    std::vector<DisplacementBC> displacement_bcs;
    std::vector<NormalDisplacementBC> normal_displacement_bcs;
    std::vector<PressureBC> pressure_bcs;
    double reference_temperature = 293.15;
};

struct FEConfig
{
    // Serial mesh kept for loading; ParMesh is always created from it.
    std::unique_ptr<mfem::Mesh> serial_mesh;
    std::unique_ptr<mfem::ParMesh> pmesh;

    std::unique_ptr<mfem::FiniteElementCollection> scalar_fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> scalar_fespace;

    std::unique_ptr<mfem::FiniteElementCollection> vector_fec;
    std::unique_ptr<mfem::ParFiniteElementSpace> vector_fespace;

    mfem::ParMesh &GetMesh()
    {
        return *pmesh;
    }
    [[nodiscard]] const mfem::ParMesh &GetMesh() const
    {
        return *pmesh;
    }
    mfem::ParFiniteElementSpace &GetScalarFESpace()
    {
        return *scalar_fespace;
    }
    [[nodiscard]] const mfem::ParFiniteElementSpace &GetScalarFESpace() const
    {
        return *scalar_fespace;
    }
    mfem::ParFiniteElementSpace &GetVectorFESpace()
    {
        return *vector_fespace;
    }
    [[nodiscard]] const mfem::ParFiniteElementSpace &GetVectorFESpace() const
    {
        return *vector_fespace;
    }
};

struct ProjectConfig
{
    SimulationConfig simulation;
    MaterialDatabase materials;

    ScalarFieldConfig electric_field;
    ScalarFieldConfig thermal_field;
    MechanicalFieldConfig mechanical_field;

    FEConfig fe;

    bool HasElectricField() const
    {
        return electric_field.enabled;
    }
    bool HasThermalField() const
    {
        return thermal_field.enabled;
    }
    bool HasMechanicalField() const
    {
        return mechanical_field.enabled;
    }
    bool NeedsPicardIteration() const;
};
}  // namespace fem::frontend

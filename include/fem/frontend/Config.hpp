#pragma once

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
    std::string output_dir = "results";
    std::string log_level = "info";
    std::string solver = "pcg";
    std::string assembly_mode = "serial";
    int picard_max_iterations = 30;
    double picard_tolerance = 1.0e-6;
    double picard_relaxation = 1.0;
    bool transient_enabled = false;
    double transient_t_start = 0.0;
    double transient_t_end = 1.0;
    double transient_dt = 1.0;
};

// 材料数据结构
struct MaterialDatabase
{
    // 建立域到材料的映射，一个域只能有一个材料，但一个材料可以对应多个域
    std::unordered_map<std::string, std::vector<int>> domain_to_material;

    // 建立材料到属性的映射，每个材料有多个属性，每个属性是一个表达式字符串
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
        material_to_properties;
};

struct ScalarFieldConfig
{
    // 标量边界条件
    struct BoundaryCondition
    {
        // 是否是Dirichlet边界条件.
        // 如果是,则q表示Dirichlet值;如果不是,则l表示Robin边界条件的系数l, q表示Robin边界条件的系数q
        bool is_dirichlet = true;
        double l = 0.0;
        double q = 0.0;
    };

    std::string type;
    std::string source_default = "0.0";
    std::unordered_map<int, std::string> domain_to_source;
    std::unordered_map<int, BoundaryCondition> boundary_to_bc;
};

struct MechanicalFieldConfig
{
    // 向量边界条件
    struct BoundaryCondition
    {
        enum class BCType
        {
            FREE,
            DIRICHLET_FULL,
            DIRICHLET_NORMAL,
            NEUMANN_PRESSURE,
        } type = BCType::FREE;  // 只需要实现这四种,不需要牵引力边界条件

        double normal_dirichlet_displacement_value = 0.0;
        std::vector<double> dirichlet_displacement_value{0.0, 0.0, 0.0};
        double pressure_value = 0.0;
    };

    std::string type;
    std::vector<double> body_force_default{0.0, 0.0, 0.0};
    std::unordered_map<int, std::vector<double>> domain_to_body_force;
    std::unordered_map<int, BoundaryCondition> boundary_to_bc;
    double reference_temperature = 293.15;  // 用于热-力耦合问题的参考温度，单位K
};

struct FEConfig
{
    mfem::Mesh *mesh = nullptr;
    mfem::FiniteElementCollection *fec = nullptr;
    mfem::FiniteElementSpace *space = nullptr;

};

struct ProjectConfig
{
    SimulationConfig simulation;
    MaterialDatabase materials;

    ScalarFieldConfig electric_field;
    ScalarFieldConfig thermal_field;
    MechanicalFieldConfig mechanical_field;

    FEConfig fe_config;
};
}  // namespace fem::frontend

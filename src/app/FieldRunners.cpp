#include "fem/app/FieldRunners.hpp"

namespace fem::app
{
int RunScalarField(const fem::frontend::ProjectConfig &config,
                   const fem::frontend::FieldConfig &field,
                   const std::vector<std::string> &coupled_variable_names,
                   fem::material::MaterialDatabase &materials);

int RunMechanicalField(const fem::frontend::ProjectConfig &config,
                       const fem::frontend::FieldConfig &field,
                       const std::vector<std::string> &coupled_variable_names,
                       fem::material::MaterialDatabase &materials);

FieldRunnerMap CreateFieldRunners()
{
    return {
        {"electrostatic", RunScalarField},
        {"thermal", RunScalarField},
        {"mechanical", RunMechanicalField},
    };
}
}  // namespace fem::app

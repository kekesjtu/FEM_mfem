#include "fem/app/FieldRunnerCommon.hpp"

#include "fem/app/FieldRunners.hpp"

namespace fem::app
{
int RunMechanicalField(const fem::frontend::ProjectConfig &config,
                       const fem::frontend::FieldConfig &field,
                       const std::vector<std::string> &coupled_variable_names,
                       fem::material::MaterialDatabase &materials)
{
    return detail::RunMechanicalFieldInternal(config, field, coupled_variable_names, materials,
                                              nullptr);
}
}  // namespace fem::app

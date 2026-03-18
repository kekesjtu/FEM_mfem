#pragma once

#include "mfem.hpp"

#include <string>

namespace fem::post
{
class SolutionTextExporter
{
  public:
    static void ExportScalarNodalTxt(const std::string &txt_path, const std::string &field_name,
                                     mfem::Mesh &mesh, const mfem::GridFunction &scalar_field,
                                     const std::string &value_name);

    static void ExportVectorNodalTxt(const std::string &txt_path, const std::string &field_name,
                                     mfem::Mesh &mesh, const mfem::GridFunction &vector_field,
                                     const std::string &prefix_name);
};
}  // namespace fem::post

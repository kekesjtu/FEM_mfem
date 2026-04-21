#pragma once

#include "mfem.hpp"

#include <string>
#include <vector>

namespace fem::output
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

    /// Export transient multi-timestep data in COMSOL-compatible wide format.
    /// Each row: x y z V@t0 T@t0 disp@t0 V@t1 T@t1 disp@t1 ...
    static void ExportTransientNodalTxt(const std::string &txt_path, mfem::Mesh &mesh,
                                        mfem::FiniteElementSpace &scalar_fespace,
                                        mfem::FiniteElementSpace &vector_fespace,
                                        const std::vector<double> &times,
                                        const std::vector<mfem::Vector> &voltages,
                                        const std::vector<mfem::Vector> &temperatures,
                                        const std::vector<mfem::Vector> &displacements);
};
}  // namespace fem::output

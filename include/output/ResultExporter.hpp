#pragma once

#include "mfem.hpp"

#include <string>
#include <vector>

namespace fem::output
{

class ResultExporter
{
  public:
    ResultExporter(const std::string &output_dir, mfem::ParMesh &mesh);

    void ExportScalar(const std::string &collection_name, const std::string &display_name,
                      mfem::ParGridFunction &gf, int cycle = 0, double time = 0.0);

    void ExportVector(const std::string &collection_name, const std::string &display_name,
                      mfem::ParGridFunction &gf, const std::string &txt_prefix = "", int cycle = 0,
                      double time = 0.0);

    void ExportTransient(mfem::ParFiniteElementSpace &scalar_fespace,
                         mfem::ParFiniteElementSpace &vector_fespace,
                         const std::vector<double> &times,
                         const std::vector<mfem::Vector> &voltages,
                         const std::vector<mfem::Vector> &temperatures,
                         const std::vector<mfem::Vector> &displacements);

  private:
    std::string output_dir_;
    mfem::ParMesh &mesh_;
};

}  // namespace fem::output

#pragma once

#include "mfem.hpp"

#include <string>
#include <vector>

namespace fem::output
{

/// Unified result exporter: combines ParaView (.vtu) + text (.txt) output.
class ResultExporter
{
  public:
    ResultExporter(const std::string &output_dir, mfem::Mesh &mesh);

    /// Export a scalar field (ParaView + nodal text).
    /// @param collection_name  ParaView collection name (e.g. "thermal")
    /// @param display_name     Field/value name for ParaView and text file (e.g. "temperature")
    void ExportScalar(const std::string &collection_name, const std::string &display_name,
                      mfem::GridFunction &gf, int cycle = 0, double time = 0.0);

    /// Export a vector field (ParaView + nodal text).
    /// @param collection_name  ParaView collection name (e.g. "mechanical")
    /// @param display_name     Field name for ParaView and text file (e.g. "displacement")
    /// @param txt_prefix       Prefix for text column headers (e.g. "u" → ux, uy, uz)
    void ExportVector(const std::string &collection_name, const std::string &display_name,
                      mfem::GridFunction &gf, const std::string &txt_prefix = "", int cycle = 0,
                      double time = 0.0);

    /// Export transient snapshots: ParaView time-series (+ txt in serial mode).
    void ExportTransient(mfem::FiniteElementSpace &scalar_fespace,
                         mfem::FiniteElementSpace &vector_fespace, const std::vector<double> &times,
                         const std::vector<mfem::Vector> &voltages,
                         const std::vector<mfem::Vector> &temperatures,
                         const std::vector<mfem::Vector> &displacements);

  private:
    std::string output_dir_;
    mfem::Mesh &mesh_;
};

}  // namespace fem::output

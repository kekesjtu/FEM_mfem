#pragma once

#include "mfem.hpp"

#include <string>

namespace fem::output
{
class ParaviewExporter
{
  public:
    explicit ParaviewExporter(std::string output_dir);

    void Save(const std::string &collection_name, const std::string &field_name, mfem::Mesh &mesh,
              mfem::GridFunction &field, int cycle, double time) const;

    void SaveTransient(const std::string &collection_name, const std::string &field_name,
                       mfem::Mesh &mesh, mfem::FiniteElementSpace &fespace,
                       const std::vector<mfem::Vector> &snapshots,
                       const std::vector<double> &times) const;

  private:
    std::string output_dir_;
};
}  // namespace fem::output

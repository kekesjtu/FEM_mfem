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

  private:
    std::string output_dir_;
};

/// Persistent ParaView collection used to write and release transient snapshots one at a time.
class ParaviewTimeSeries
{
  public:
    ParaviewTimeSeries(const std::string &output_dir, const std::string &collection_name,
                       const std::string &field_name, mfem::Mesh &mesh,
                       mfem::FiniteElementSpace &fespace);

    void SaveSnapshot(const mfem::Vector &snapshot, int cycle, double time);

  private:
    mfem::GridFunction field_;
    mfem::ParaViewDataCollection collection_;
};
}  // namespace fem::output

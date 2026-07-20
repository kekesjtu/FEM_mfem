#pragma once

#include "mfem.hpp"

#include <cstddef>
#include <memory>
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
};

/// Accumulates only the root-rank nodal values needed by the existing
/// COMSOL-compatible transient wide-table format. Serial GridFunctions can be
/// consumed and released one snapshot at a time.
class TransientNodalTextWriter
{
  public:
    TransientNodalTextWriter(std::string txt_path, mfem::Mesh &mesh,
                             const std::vector<double> &times, bool has_voltage,
                             bool has_temperature, bool has_displacement);
    ~TransientNodalTextWriter();

    TransientNodalTextWriter(const TransientNodalTextWriter &) = delete;
    TransientNodalTextWriter &operator=(const TransientNodalTextWriter &) = delete;
    TransientNodalTextWriter(TransientNodalTextWriter &&) noexcept;
    TransientNodalTextWriter &operator=(TransientNodalTextWriter &&) noexcept;

    void AddVoltage(std::size_t snapshot_index, const mfem::GridFunction &field);
    void AddTemperature(std::size_t snapshot_index, const mfem::GridFunction &field);
    void AddDisplacement(std::size_t snapshot_index, const mfem::GridFunction &field);
    void Write() const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}  // namespace fem::output

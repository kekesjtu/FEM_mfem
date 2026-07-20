#include "output/ParaviewExporter.hpp"

#include <filesystem>

namespace fem::output
{
ParaviewExporter::ParaviewExporter(std::string output_dir) : output_dir_(std::move(output_dir))
{
}

void ParaviewExporter::Save(const std::string &collection_name, const std::string &field_name,
                            mfem::Mesh &mesh, mfem::GridFunction &field, int cycle,
                            double time) const
{
    std::filesystem::create_directories(output_dir_);

    mfem::ParaViewDataCollection collection(collection_name, &mesh);
    collection.SetPrefixPath(output_dir_);
    collection.SetLevelsOfDetail(1);
    collection.SetCycle(cycle);
    collection.SetTime(time);
    collection.RegisterField(field_name, &field);
    collection.SetHighOrderOutput(true);
    collection.SetDataFormat(mfem::VTKFormat::BINARY);
    collection.Save();
}

ParaviewTimeSeries::ParaviewTimeSeries(const std::string &output_dir,
                                       const std::string &collection_name,
                                       const std::string &field_name, mfem::Mesh &mesh,
                                       mfem::FiniteElementSpace &fespace)
    : field_(&fespace), collection_(collection_name, &mesh)
{
    std::filesystem::create_directories(output_dir);
    collection_.SetPrefixPath(output_dir);
    collection_.SetLevelsOfDetail(1);
    collection_.SetHighOrderOutput(true);
    collection_.SetDataFormat(mfem::VTKFormat::BINARY);
    collection_.RegisterField(field_name, &field_);
}

void ParaviewTimeSeries::SaveSnapshot(const mfem::Vector &snapshot, int cycle, double time)
{
    field_ = snapshot;
    collection_.SetCycle(cycle);
    collection_.SetTime(time);
    collection_.Save();
}
}  // namespace fem::output

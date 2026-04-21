#include "fem/output/ParaviewExporter.hpp"

#include <filesystem>
#include <stdexcept>

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

void ParaviewExporter::SaveTransient(const std::string &collection_name,
                                     const std::string &field_name, mfem::Mesh &mesh,
                                     mfem::FiniteElementSpace &fespace,
                                     const std::vector<mfem::Vector> &snapshots,
                                     const std::vector<double> &times) const
{
    if (snapshots.size() != times.size())
    {
        throw std::runtime_error(
            "ParaviewExporter::SaveTransient requires snapshots/time size match.");
    }
    if (snapshots.empty())
    {
        return;
    }

    std::filesystem::create_directories(output_dir_);

    mfem::GridFunction gf(&fespace);
    mfem::ParaViewDataCollection collection(collection_name, &mesh);
    collection.SetPrefixPath(output_dir_);
    collection.SetLevelsOfDetail(1);
    collection.SetHighOrderOutput(true);
    collection.SetDataFormat(mfem::VTKFormat::BINARY);
    collection.RegisterField(field_name, &gf);

    for (size_t i = 0; i < snapshots.size(); ++i)
    {
        gf = snapshots[i];
        collection.SetCycle(static_cast<int>(i));
        collection.SetTime(times[i]);
        collection.Save();
    }
}
}  // namespace fem::output

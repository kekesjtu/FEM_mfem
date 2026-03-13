#include "fem/io/ParaviewExporter.hpp"

#include <filesystem>

namespace fem::io
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
    collection.Save();
}
}  // namespace fem::io

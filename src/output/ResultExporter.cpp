#include "fem/output/ResultExporter.hpp"

#include "fem/output/ParaviewExporter.hpp"
#include "fem/output/SolutionTextExporter.hpp"

#include <filesystem>

namespace fem::io
{

namespace
{
bool IsParallelMesh([[maybe_unused]] mfem::Mesh &mesh)
{
#ifdef MFEM_USE_MPI
    return dynamic_cast<mfem::ParMesh *>(&mesh) != nullptr;
#else
    return false;
#endif
}
}  // namespace

ResultExporter::ResultExporter(const std::string &output_dir, mfem::Mesh &mesh)
    : output_dir_(output_dir), mesh_(mesh)
{
    std::filesystem::create_directories(output_dir_);
}

void ResultExporter::ExportScalar(const std::string &collection_name,
                                  const std::string &display_name, mfem::GridFunction &gf,
                                  int cycle, double time)
{
    ParaviewExporter pv(output_dir_);
    pv.Save(collection_name, display_name, mesh_, gf, cycle, time);

    // txt export only in serial — parallel results use ParaView for comparison
    if (!IsParallelMesh(mesh_))
    {
        post::SolutionTextExporter::ExportScalarNodalTxt(output_dir_ + "/" + display_name + ".txt",
                                                         collection_name, mesh_, gf, display_name);
    }
}

void ResultExporter::ExportVector(const std::string &collection_name,
                                  const std::string &display_name, mfem::GridFunction &gf,
                                  const std::string &txt_prefix, int cycle, double time)
{
    ParaviewExporter pv(output_dir_);
    pv.Save(collection_name, display_name, mesh_, gf, cycle, time);

    if (!IsParallelMesh(mesh_))
    {
        const std::string &prefix = txt_prefix.empty() ? display_name : txt_prefix;
        post::SolutionTextExporter::ExportVectorNodalTxt(output_dir_ + "/" + display_name + ".txt",
                                                         collection_name, mesh_, gf, prefix);
    }
}

}  // namespace fem::io

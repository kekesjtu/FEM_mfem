#include "fem/output/ResultExporter.hpp"

#include "fem/output/ParaviewExporter.hpp"
#include "fem/output/SolutionTextExporter.hpp"

#include <filesystem>

namespace fem::output
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
        SolutionTextExporter::ExportScalarNodalTxt(output_dir_ + "/" + display_name + ".txt",
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
        SolutionTextExporter::ExportVectorNodalTxt(output_dir_ + "/" + display_name + ".txt",
                                                   collection_name, mesh_, gf, prefix);
    }
}

void ResultExporter::ExportTransient(mfem::FiniteElementSpace &scalar_fespace,
                                     mfem::FiniteElementSpace &vector_fespace,
                                     const std::vector<double> &times,
                                     const std::vector<mfem::Vector> &voltages,
                                     const std::vector<mfem::Vector> &temperatures,
                                     const std::vector<mfem::Vector> &displacements)
{
    // ParaView time-series export
    ParaviewExporter pv(output_dir_);
    if (!voltages.empty() && voltages.front().Size() > 0)
    {
        pv.SaveTransient("electrostatic", "voltage", mesh_, scalar_fespace, voltages, times);
    }
    if (!temperatures.empty() && temperatures.front().Size() > 0)
    {
        pv.SaveTransient("thermal", "temperature", mesh_, scalar_fespace, temperatures, times);
    }
    if (!displacements.empty() && displacements.front().Size() > 0)
    {
        pv.SaveTransient("mechanical", "displacement", mesh_, vector_fespace, displacements, times);
    }

    // txt export only in serial mode
    if (!IsParallelMesh(mesh_))
    {
        SolutionTextExporter::ExportTransientNodalTxt(output_dir_ + "/transient.txt", mesh_,
                                                      scalar_fespace, vector_fespace, times,
                                                      voltages, temperatures, displacements);
    }
}

}  // namespace fem::output

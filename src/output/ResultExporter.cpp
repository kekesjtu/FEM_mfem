#include "fem/output/ResultExporter.hpp"

#include "fem/output/ParaviewExporter.hpp"
#include "fem/output/SolutionTextExporter.hpp"

#include <mpi.h>
#include <filesystem>

namespace fem::output
{

ResultExporter::ResultExporter(const std::string &output_dir, mfem::ParMesh &mesh)
    : output_dir_(output_dir), mesh_(mesh)
{
    int rank, np;
    MPI_Comm_rank(mesh.GetComm(), &rank);
    MPI_Comm_size(mesh.GetComm(), &np);

    if (rank == 0)
        std::filesystem::remove_all(output_dir_);
    if (np > 1)
        MPI_Barrier(mesh.GetComm());
    std::filesystem::create_directories(output_dir_);
}

void ResultExporter::ExportScalar(const std::string &collection_name,
                                  const std::string &display_name, mfem::ParGridFunction &gf,
                                  int cycle, double time)
{
    ParaviewExporter pv(output_dir_);
    pv.Save(collection_name, display_name, mesh_, gf, cycle, time);

    int rank;
    MPI_Comm_rank(mesh_.GetComm(), &rank);
    auto serial_mesh = mesh_.GetSerialMesh(0);
    auto serial_gf = gf.GetSerialGridFunction(0, serial_mesh);
    if (rank == 0)
    {
        SolutionTextExporter::ExportScalarNodalTxt(output_dir_ + "/" + display_name + ".txt",
                                                   collection_name, serial_mesh, serial_gf,
                                                   display_name);
    }
}

void ResultExporter::ExportVector(const std::string &collection_name,
                                  const std::string &display_name, mfem::ParGridFunction &gf,
                                  const std::string &txt_prefix, int cycle, double time)
{
    ParaviewExporter pv(output_dir_);
    pv.Save(collection_name, display_name, mesh_, gf, cycle, time);

    int rank;
    MPI_Comm_rank(mesh_.GetComm(), &rank);
    auto serial_mesh = mesh_.GetSerialMesh(0);
    auto serial_gf = gf.GetSerialGridFunction(0, serial_mesh);
    if (rank == 0)
    {
        const std::string &prefix = txt_prefix.empty() ? display_name : txt_prefix;
        SolutionTextExporter::ExportVectorNodalTxt(output_dir_ + "/" + display_name + ".txt",
                                                   collection_name, serial_mesh, serial_gf, prefix);
    }
}

void ResultExporter::ExportTransient(mfem::ParFiniteElementSpace &scalar_fespace,
                                     mfem::ParFiniteElementSpace &vector_fespace,
                                     const std::vector<double> &times,
                                     const std::vector<mfem::Vector> &voltages,
                                     const std::vector<mfem::Vector> &temperatures,
                                     const std::vector<mfem::Vector> &displacements)
{
    ParaviewExporter pv(output_dir_);
    if (!voltages.empty() && voltages.front().Size() > 0)
        pv.SaveTransient("electrostatic", "voltage", mesh_, scalar_fespace, voltages, times);
    if (!temperatures.empty() && temperatures.front().Size() > 0)
        pv.SaveTransient("thermal", "temperature", mesh_, scalar_fespace, temperatures, times);
    if (!displacements.empty() && displacements.front().Size() > 0)
        pv.SaveTransient("mechanical", "displacement", mesh_, vector_fespace, displacements, times);

    int rank;
    MPI_Comm_rank(mesh_.GetComm(), &rank);
    auto serial_mesh = mesh_.GetSerialMesh(0);

    // Build serial FE spaces for txt export on rank 0
    const int dim = mesh_.Dimension();
    const int order = scalar_fespace.FEColl()->GetOrder();
    mfem::H1_FECollection serial_scalar_fec(order, dim);
    mfem::H1_FECollection serial_vector_fec(order, dim);
    mfem::FiniteElementSpace serial_scalar_fes(&serial_mesh, &serial_scalar_fec);
    mfem::FiniteElementSpace serial_vector_fes(&serial_mesh, &serial_vector_fec, dim);

    // Gather each snapshot to serial GridFunction
    auto gather = [&](const std::vector<mfem::Vector> &par_data,
                      mfem::ParFiniteElementSpace &pfes) -> std::vector<mfem::Vector>
    {
        std::vector<mfem::Vector> result;
        mfem::ParGridFunction pgf(&pfes);
        for (const auto &v : par_data)
        {
            pgf = v;
            auto sgf = pgf.GetSerialGridFunction(0, serial_mesh);
            if (rank == 0)
                result.push_back(sgf);
            else
                result.emplace_back();
        }
        return result;
    };

    std::vector<mfem::Vector> sv, st, sd;
    if (!voltages.empty() && voltages.front().Size() > 0)
        sv = gather(voltages, scalar_fespace);
    if (!temperatures.empty() && temperatures.front().Size() > 0)
        st = gather(temperatures, scalar_fespace);
    if (!displacements.empty() && displacements.front().Size() > 0)
        sd = gather(displacements, vector_fespace);

    if (rank == 0)
    {
        SolutionTextExporter::ExportTransientNodalTxt(output_dir_ + "/transient.txt", serial_mesh,
                                                      serial_scalar_fes, serial_vector_fes, times,
                                                      sv, st, sd);
    }
}

}  // namespace fem::output

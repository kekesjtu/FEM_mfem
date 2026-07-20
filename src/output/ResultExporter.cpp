#include "output/ResultExporter.hpp"

#include "output/ParaviewExporter.hpp"
#include "output/SolutionTextExporter.hpp"

#include <mpi.h>
#include <filesystem>
#include <memory>
#include <stdexcept>

namespace
{
bool HasFieldSnapshots(const char *field_name, std::size_t expected_count,
                       const std::vector<mfem::Vector> &snapshots)
{
    if (snapshots.empty())
        return false;
    if (snapshots.size() != expected_count)
    {
        throw std::runtime_error(std::string("Transient ") + field_name +
                                 " snapshot/time size mismatch");
    }

    const int field_size = snapshots.front().Size();
    for (const auto &snapshot : snapshots)
    {
        if (snapshot.Size() != field_size)
        {
            throw std::runtime_error(std::string("Transient ") + field_name +
                                     " snapshots have inconsistent sizes");
        }
    }
    return field_size > 0;
}
}  // namespace

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

void ResultExporter::ExportTransient(mfem::ParFiniteElementSpace *scalar_fespace,
                                     mfem::ParFiniteElementSpace *vector_fespace,
                                     const std::vector<double> &times,
                                     std::vector<mfem::Vector> &voltages,
                                     std::vector<mfem::Vector> &temperatures,
                                     std::vector<mfem::Vector> &displacements)
{
    const bool has_voltage = HasFieldSnapshots("voltage", times.size(), voltages);
    const bool has_temperature = HasFieldSnapshots("temperature", times.size(), temperatures);
    const bool has_displacement = HasFieldSnapshots("displacement", times.size(), displacements);

    if ((has_voltage || has_temperature) && !scalar_fespace)
        throw std::runtime_error("Transient scalar output requires a scalar finite element space");
    if (has_displacement && !vector_fespace)
        throw std::runtime_error(
            "Transient displacement output requires a vector finite element space");

    std::unique_ptr<ParaviewTimeSeries> voltage_series;
    std::unique_ptr<ParaviewTimeSeries> temperature_series;
    std::unique_ptr<ParaviewTimeSeries> displacement_series;
    if (has_voltage)
        voltage_series = std::make_unique<ParaviewTimeSeries>(output_dir_, "electrostatic",
                                                              "voltage", mesh_, *scalar_fespace);
    if (has_temperature)
        temperature_series = std::make_unique<ParaviewTimeSeries>(
            output_dir_, "thermal", "temperature", mesh_, *scalar_fespace);
    if (has_displacement)
        displacement_series = std::make_unique<ParaviewTimeSeries>(
            output_dir_, "mechanical", "displacement", mesh_, *vector_fespace);

    int rank;
    MPI_Comm_rank(mesh_.GetComm(), &rank);
    auto serial_mesh = mesh_.GetSerialMesh(0);

    std::unique_ptr<mfem::FiniteElementCollection> serial_scalar_fec;
    std::unique_ptr<mfem::FiniteElementCollection> serial_vector_fec;
    std::unique_ptr<mfem::FiniteElementSpace> serial_scalar_fes;
    std::unique_ptr<mfem::FiniteElementSpace> serial_vector_fes;
    std::unique_ptr<mfem::ParGridFunction> scalar_field;
    std::unique_ptr<mfem::ParGridFunction> vector_field;

    if (has_voltage || has_temperature)
    {
        serial_scalar_fec.reset(
            scalar_fespace->FEColl()->Clone(scalar_fespace->FEColl()->GetOrder()));
        serial_scalar_fes = std::make_unique<mfem::FiniteElementSpace>(
            &serial_mesh, serial_scalar_fec.get(), scalar_fespace->GetVDim(),
            scalar_fespace->GetOrdering());
        scalar_field = std::make_unique<mfem::ParGridFunction>(scalar_fespace);
    }
    if (has_displacement)
    {
        serial_vector_fec.reset(
            vector_fespace->FEColl()->Clone(vector_fespace->FEColl()->GetOrder()));
        serial_vector_fes = std::make_unique<mfem::FiniteElementSpace>(
            &serial_mesh, serial_vector_fec.get(), vector_fespace->GetVDim(),
            vector_fespace->GetOrdering());
        vector_field = std::make_unique<mfem::ParGridFunction>(vector_fespace);
    }

    std::unique_ptr<TransientNodalTextWriter> text_writer;
    if (rank == 0)
    {
        text_writer = std::make_unique<TransientNodalTextWriter>(output_dir_ + "/transient.txt",
                                                                 serial_mesh, times, has_voltage,
                                                                 has_temperature, has_displacement);
    }

    for (std::size_t i = 0; i < times.size(); ++i)
    {
        const int cycle = static_cast<int>(i);
        if (has_voltage)
        {
            voltage_series->SaveSnapshot(voltages[i], cycle, times[i]);
            *scalar_field = voltages[i];
            auto serial_field = scalar_field->GetSerialGridFunction(0, *serial_scalar_fes);
            if (rank == 0)
                text_writer->AddVoltage(i, serial_field);
            voltages[i].Destroy();
        }
        if (has_temperature)
        {
            temperature_series->SaveSnapshot(temperatures[i], cycle, times[i]);
            *scalar_field = temperatures[i];
            auto serial_field = scalar_field->GetSerialGridFunction(0, *serial_scalar_fes);
            if (rank == 0)
                text_writer->AddTemperature(i, serial_field);
            temperatures[i].Destroy();
        }
        if (has_displacement)
        {
            displacement_series->SaveSnapshot(displacements[i], cycle, times[i]);
            *vector_field = displacements[i];
            auto serial_field = vector_field->GetSerialGridFunction(0, *serial_vector_fes);
            if (rank == 0)
                text_writer->AddDisplacement(i, serial_field);
            displacements[i].Destroy();
        }
    }

    voltages.clear();
    temperatures.clear();
    displacements.clear();

    if (rank == 0)
        text_writer->Write();
}

}  // namespace fem::output

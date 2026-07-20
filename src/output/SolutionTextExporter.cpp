#include "output/SolutionTextExporter.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fem::output
{
namespace
{

struct VertexSampleLocation
{
    int element_id = -1;
    int local_vertex_id = -1;
};

struct VertexCoord
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

std::vector<VertexSampleLocation> BuildVertexSampleLocations(mfem::Mesh &mesh)
{
    std::vector<VertexSampleLocation> locations(mesh.GetNV());
    mfem::Array<int> vertices;
    for (int e = 0; e < mesh.GetNE(); ++e)
    {
        mesh.GetElementVertices(e, vertices);
        for (int lv = 0; lv < vertices.Size(); ++lv)
        {
            const int v = vertices[lv];
            if (v < 0 || v >= static_cast<int>(locations.size()))
                continue;
            if (locations[v].element_id < 0)
            {
                locations[v].element_id = e;
                locations[v].local_vertex_id = lv;
            }
        }
    }
    return locations;
}

mfem::IntegrationPoint VertexIntegrationPoint(const mfem::Mesh &mesh, int element_id,
                                              int local_vertex_id)
{
    const int geom = mesh.GetElementBaseGeometry(element_id);
    const mfem::IntegrationRule *ir = mfem::Geometries.GetVertices(geom);
    return ir->IntPoint(local_vertex_id);
}

}  // namespace

void SolutionTextExporter::ExportScalarNodalTxt(const std::string &txt_path,
                                                const std::string &field_name, mfem::Mesh &mesh,
                                                const mfem::GridFunction &scalar_field,
                                                const std::string &value_name)
{
    const int nv = mesh.GetNV();
    const int dim = mesh.Dimension();
    const auto locations = BuildVertexSampleLocations(mesh);

    std::filesystem::create_directories(std::filesystem::path(txt_path).parent_path());
    std::ofstream out(txt_path);
    if (!out)
        throw std::runtime_error("Failed to open output txt file: " + txt_path);
    out << std::setprecision(17);

    out << "% Model:              FEM_mfem\n"
        << "% Field:              " << field_name << "\n"
        << "% Dimension:          " << dim << "\n"
        << "% Nodes:              " << nv << "\n"
        << "% Expressions:        1\n"
        << "% Description:        " << value_name << "\n"
        << "% x                       y                        z                        "
        << value_name << "\n";

    mfem::Vector pt(dim);
    for (int v = 0; v < nv; ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
            continue;
        auto *tr = mesh.GetElementTransformation(loc.element_id);
        auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
        tr->SetIntPoint(&ip);
        tr->Transform(ip, pt);
        double val = scalar_field.GetValue(*tr, ip);
        out << pt(0) << " " << pt(1) << " " << (dim > 2 ? pt(2) : 0.0) << " " << val << "\n";
    }
}

void SolutionTextExporter::ExportVectorNodalTxt(const std::string &txt_path,
                                                const std::string &field_name, mfem::Mesh &mesh,
                                                const mfem::GridFunction &vector_field,
                                                const std::string &prefix_name)
{
    const int nv = mesh.GetNV();
    const int dim = mesh.Dimension();
    const auto locations = BuildVertexSampleLocations(mesh);

    std::filesystem::create_directories(std::filesystem::path(txt_path).parent_path());
    std::ofstream out(txt_path);
    if (!out)
        throw std::runtime_error("Failed to open output txt file: " + txt_path);
    out << std::setprecision(17);

    out << "% Model:              FEM_mfem\n"
        << "% Field:              " << field_name << "\n"
        << "% Dimension:          " << dim << "\n"
        << "% Nodes:              " << nv << "\n"
        << "% Expressions:        " << (dim + 1) << "\n"
        << "% Description:        " << prefix_name << " components and magnitude\n"
        << "% x                       y                        z                        "
        << prefix_name << "x  " << prefix_name << "y  " << prefix_name << "z  " << prefix_name
        << "_magnitude\n";

    mfem::Vector vec(dim);
    mfem::Vector pt(dim);
    for (int v = 0; v < nv; ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
            continue;
        auto *tr = mesh.GetElementTransformation(loc.element_id);
        auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
        tr->SetIntPoint(&ip);
        tr->Transform(ip, pt);
        vector_field.GetVectorValue(*tr, ip, vec);
        double mag2 = 0.0;
        for (int d = 0; d < dim; ++d)
            mag2 += vec(d) * vec(d);
        out << pt(0) << " " << pt(1) << " " << (dim > 2 ? pt(2) : 0.0);
        for (int d = 0; d < dim; ++d)
            out << " " << vec(d);
        out << " " << std::sqrt(mag2) << "\n";
    }
}

struct TransientNodalTextWriter::Impl
{
    Impl(std::string path, mfem::Mesh &serial_mesh, const std::vector<double> &snapshot_times,
         bool voltage_enabled, bool temperature_enabled, bool displacement_enabled)
        : txt_path(std::move(path)),
          mesh(serial_mesh),
          times(snapshot_times),
          has_voltage(voltage_enabled),
          has_temperature(temperature_enabled),
          has_displacement(displacement_enabled),
          nv(mesh.GetNV()),
          dim(mesh.Dimension()),
          locations(BuildVertexSampleLocations(mesh)),
          coords(static_cast<std::size_t>(nv))
    {
        int next_offset = 0;
        if (has_voltage)
            voltage_offset = next_offset++;
        if (has_temperature)
            temperature_offset = next_offset++;
        if (has_displacement)
            displacement_offset = next_offset++;
        nfields = next_offset;

        const auto snapshot_count = times.size();
        const auto vertex_count = static_cast<std::size_t>(nv);
        const auto field_count = static_cast<std::size_t>(nfields);
        if (field_count > 0 && vertex_count > 0 &&
            snapshot_count > std::numeric_limits<std::size_t>::max() / vertex_count / field_count)
        {
            throw std::length_error("Transient nodal output buffer size overflow");
        }
        data.resize(snapshot_count * vertex_count * field_count);

        mfem::Vector point(dim);
        for (int v = 0; v < nv; ++v)
        {
            const auto &loc = locations[v];
            if (loc.element_id < 0)
                continue;
            auto *tr = mesh.GetElementTransformation(loc.element_id);
            auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
            tr->SetIntPoint(&ip);
            tr->Transform(ip, point);
            coords[static_cast<std::size_t>(v)] = {point(0), point(1), dim > 2 ? point(2) : 0.0};
        }
    }

    std::size_t DataIndex(std::size_t snapshot_index, int vertex, int field_offset) const
    {
        return (snapshot_index * static_cast<std::size_t>(nv) + static_cast<std::size_t>(vertex)) *
                   static_cast<std::size_t>(nfields) +
               static_cast<std::size_t>(field_offset);
    }

    void ValidateSnapshot(std::size_t snapshot_index, int field_offset,
                          const char *field_name) const
    {
        if (snapshot_index >= times.size())
            throw std::out_of_range("Transient snapshot index is out of range");
        if (field_offset < 0)
            throw std::logic_error(std::string("Transient output does not contain field: ") +
                                   field_name);
    }

    void AddScalar(std::size_t snapshot_index, const mfem::GridFunction &field, int field_offset,
                   const char *field_name)
    {
        ValidateSnapshot(snapshot_index, field_offset, field_name);
        for (int v = 0; v < nv; ++v)
        {
            const auto &loc = locations[static_cast<std::size_t>(v)];
            if (loc.element_id < 0)
                continue;
            auto *tr = mesh.GetElementTransformation(loc.element_id);
            auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
            tr->SetIntPoint(&ip);
            data[DataIndex(snapshot_index, v, field_offset)] = field.GetValue(*tr, ip);
        }
    }

    void AddVectorMagnitude(std::size_t snapshot_index, const mfem::GridFunction &field)
    {
        ValidateSnapshot(snapshot_index, displacement_offset, "displacement");
        mfem::Vector value(dim);
        for (int v = 0; v < nv; ++v)
        {
            const auto &loc = locations[static_cast<std::size_t>(v)];
            if (loc.element_id < 0)
                continue;
            auto *tr = mesh.GetElementTransformation(loc.element_id);
            auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
            tr->SetIntPoint(&ip);
            field.GetVectorValue(*tr, ip, value);

            double magnitude_sq = 0.0;
            for (int d = 0; d < dim; ++d)
                magnitude_sq += value(d) * value(d);
            data[DataIndex(snapshot_index, v, displacement_offset)] = std::sqrt(magnitude_sq);
        }
    }

    void Write() const
    {
        std::filesystem::create_directories(std::filesystem::path(txt_path).parent_path());
        std::ofstream out(txt_path);
        if (!out)
            throw std::runtime_error("Failed to open output txt file: " + txt_path);
        out << std::setprecision(17);

        out << "% Model:              FEM_mfem\n"
            << "% Dimension:          " << dim << "\n"
            << "% Nodes:              " << nv << "\n"
            << "% Time steps:         " << times.size() << "\n"
            << "% Description:        Electric potential, Temperature, Displacement magnitude\n"
            << "% x                       y                        z                        ";
        for (std::size_t s = 0; s < times.size(); ++s)
        {
            if (has_voltage)
                out << "V @ t=" << times[s] << "  ";
            if (has_temperature)
                out << "T @ t=" << times[s] << "  ";
            if (has_displacement)
                out << "disp @ t=" << times[s] << "  ";
        }
        out << "\n";

        for (int v = 0; v < nv; ++v)
        {
            const auto &loc = locations[static_cast<std::size_t>(v)];
            if (loc.element_id < 0)
                continue;
            const auto &coord = coords[static_cast<std::size_t>(v)];
            out << coord.x << " " << coord.y << " " << coord.z;
            for (std::size_t s = 0; s < times.size(); ++s)
            {
                for (int field_offset = 0; field_offset < nfields; ++field_offset)
                    out << " " << data[DataIndex(s, v, field_offset)];
            }
            out << "\n";
        }
    }

    std::string txt_path;
    mfem::Mesh &mesh;
    std::vector<double> times;
    bool has_voltage = false;
    bool has_temperature = false;
    bool has_displacement = false;
    int nv = 0;
    int dim = 0;
    int nfields = 0;
    int voltage_offset = -1;
    int temperature_offset = -1;
    int displacement_offset = -1;
    std::vector<VertexSampleLocation> locations;
    std::vector<VertexCoord> coords;
    std::vector<double> data;
};

TransientNodalTextWriter::TransientNodalTextWriter(std::string txt_path, mfem::Mesh &mesh,
                                                   const std::vector<double> &times,
                                                   bool has_voltage, bool has_temperature,
                                                   bool has_displacement)
    : impl_(std::make_unique<Impl>(std::move(txt_path), mesh, times, has_voltage, has_temperature,
                                   has_displacement))
{
}

TransientNodalTextWriter::~TransientNodalTextWriter() = default;
TransientNodalTextWriter::TransientNodalTextWriter(TransientNodalTextWriter &&) noexcept = default;
TransientNodalTextWriter &TransientNodalTextWriter::operator=(
    TransientNodalTextWriter &&) noexcept = default;

void TransientNodalTextWriter::AddVoltage(std::size_t snapshot_index,
                                          const mfem::GridFunction &field)
{
    impl_->AddScalar(snapshot_index, field, impl_->voltage_offset, "voltage");
}

void TransientNodalTextWriter::AddTemperature(std::size_t snapshot_index,
                                              const mfem::GridFunction &field)
{
    impl_->AddScalar(snapshot_index, field, impl_->temperature_offset, "temperature");
}

void TransientNodalTextWriter::AddDisplacement(std::size_t snapshot_index,
                                               const mfem::GridFunction &field)
{
    impl_->AddVectorMagnitude(snapshot_index, field);
}

void TransientNodalTextWriter::Write() const
{
    impl_->Write();
}

}  // namespace fem::output

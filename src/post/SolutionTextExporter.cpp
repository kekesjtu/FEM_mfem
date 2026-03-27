#include "fem/post/SolutionTextExporter.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace fem::post
{
namespace
{
std::string ParallelAwareOutputPath(const std::string &path)
{
#if defined(MFEM_USE_MPI)
    if (mfem::Mpi::IsInitialized() && mfem::Mpi::WorldSize() > 1)
    {
        const int rank = mfem::Mpi::WorldRank();
        const std::string suffix = ".rank" + std::to_string(rank);
        const std::size_t dot = path.find_last_of('.');
        if (dot == std::string::npos)
        {
            return path + suffix;
        }
        return path.substr(0, dot) + suffix + path.substr(dot);
    }
#endif
    return path;
}

struct VertexSampleLocation
{
    int element_id = -1;
    int local_vertex_id = -1;
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
            {
                continue;
            }
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
    const mfem::IntegrationRule *vertices = mfem::Geometries.GetVertices(geom);
    if (!vertices || local_vertex_id < 0 || local_vertex_id >= vertices->GetNPoints())
    {
        throw std::runtime_error("Invalid local vertex id when exporting nodal txt.");
    }
    return vertices->IntPoint(local_vertex_id);
}

void WriteScalarHeader(std::ofstream &out, const std::string &field_name, int dim, int nodes,
                       const std::string &value_name)
{
    out << "% Model:              FEM_mfem\n";
    out << "% Field:              " << field_name << "\n";
    out << "% Dimension:          " << dim << "\n";
    out << "% Nodes:              " << nodes << "\n";
    out << "% Expressions:        1\n";
    out << "% Description:        " << value_name << "\n";
    out << "% x                       y                        z                        "
        << value_name << "\n";
}

void WriteVectorHeader(std::ofstream &out, const std::string &field_name, int dim, int nodes,
                       const std::string &prefix_name)
{
    out << "% Model:              FEM_mfem\n";
    out << "% Field:              " << field_name << "\n";
    out << "% Dimension:          " << dim << "\n";
    out << "% Nodes:              " << nodes << "\n";
    out << "% Expressions:        " << (dim + 1) << "\n";
    out << "% Description:        " << prefix_name << " components and magnitude\n";
    out << "% x                       y                        z                        "
        << prefix_name << "x                      " << prefix_name << "y                      "
        << prefix_name << "z                      " << prefix_name << "_magnitude\n";
}
}  // namespace

void SolutionTextExporter::ExportScalarNodalTxt(const std::string &txt_path,
                                                const std::string &field_name, mfem::Mesh &mesh,
                                                const mfem::GridFunction &scalar_field,
                                                const std::string &value_name)
{
    const std::string resolved_path = ParallelAwareOutputPath(txt_path);
    std::filesystem::create_directories(std::filesystem::path(resolved_path).parent_path());

    std::ofstream out(resolved_path);
    if (!out)
    {
        throw std::runtime_error("Failed to open output txt file: " + resolved_path);
    }
    out << std::setprecision(17);

    const int dim = mesh.Dimension();
    WriteScalarHeader(out, field_name, dim, mesh.GetNV(), value_name);

    const auto locations = BuildVertexSampleLocations(mesh);

    for (int v = 0; v < mesh.GetNV(); ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
        {
            continue;
        }

        auto *tr = mesh.GetElementTransformation(loc.element_id);
        const mfem::IntegrationPoint ip =
            VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
        tr->SetIntPoint(&ip);

        mfem::Vector point;
        tr->Transform(ip, point);
        const double x = point.Size() > 0 ? point(0) : 0.0;
        const double y = point.Size() > 1 ? point(1) : 0.0;
        const double z = point.Size() > 2 ? point(2) : 0.0;
        const double value = scalar_field.GetValue(*tr, ip);

        out << x << " " << y << " " << z << " " << value << "\n";
    }
}

void SolutionTextExporter::ExportVectorNodalTxt(const std::string &txt_path,
                                                const std::string &field_name, mfem::Mesh &mesh,
                                                const mfem::GridFunction &vector_field,
                                                const std::string &prefix_name)
{
    const std::string resolved_path = ParallelAwareOutputPath(txt_path);
    std::filesystem::create_directories(std::filesystem::path(resolved_path).parent_path());

    std::ofstream out(resolved_path);
    if (!out)
    {
        throw std::runtime_error("Failed to open output txt file: " + resolved_path);
    }
    out << std::setprecision(17);

    const int dim = mesh.Dimension();
    WriteVectorHeader(out, field_name, dim, mesh.GetNV(), prefix_name);

    const auto locations = BuildVertexSampleLocations(mesh);

    mfem::Vector vec(dim);
    for (int v = 0; v < mesh.GetNV(); ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
        {
            continue;
        }

        auto *tr = mesh.GetElementTransformation(loc.element_id);
        const mfem::IntegrationPoint ip =
            VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
        tr->SetIntPoint(&ip);

        mfem::Vector point;
        tr->Transform(ip, point);
        const double x = point.Size() > 0 ? point(0) : 0.0;
        const double y = point.Size() > 1 ? point(1) : 0.0;
        const double z = point.Size() > 2 ? point(2) : 0.0;

        vector_field.GetVectorValue(*tr, ip, vec);
        const double ux = vec.Size() > 0 ? vec(0) : 0.0;
        const double uy = vec.Size() > 1 ? vec(1) : 0.0;
        const double uz = vec.Size() > 2 ? vec(2) : 0.0;
        const double umag = std::sqrt(ux * ux + uy * uy + uz * uz);

        out << x << " " << y << " " << z << " " << ux << " " << uy << " " << uz << " " << umag
            << "\n";
    }
}
}  // namespace fem::post

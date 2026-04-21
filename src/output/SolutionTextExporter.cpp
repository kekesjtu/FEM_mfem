#include "fem/output/SolutionTextExporter.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
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

void SolutionTextExporter::ExportTransientNodalTxt(const std::string &txt_path, mfem::Mesh &mesh,
                                                   mfem::FiniteElementSpace &scalar_fespace,
                                                   mfem::FiniteElementSpace &vector_fespace,
                                                   const std::vector<double> &times,
                                                   const std::vector<mfem::Vector> &voltages,
                                                   const std::vector<mfem::Vector> &temperatures,
                                                   const std::vector<mfem::Vector> &displacements)
{
    const int nsnaps = static_cast<int>(times.size());
    const int nv = mesh.GetNV();
    const int dim = mesh.Dimension();
    const bool has_v = !voltages.empty() && voltages[0].Size() > 0;
    const bool has_t = !temperatures.empty() && temperatures[0].Size() > 0;
    const bool has_d = !displacements.empty() && displacements[0].Size() > 0;
    const int fields_per_snap = (has_v ? 1 : 0) + (has_t ? 1 : 0) + (has_d ? 1 : 0);

    const auto locations = BuildVertexSampleLocations(mesh);

    // Create reusable temporary GridFunctions (one per field type)
    std::unique_ptr<mfem::GridFunction> v_gf, t_gf, d_gf;
    if (has_v)
        v_gf = std::make_unique<mfem::GridFunction>(&scalar_fespace);
    if (has_t)
        t_gf = std::make_unique<mfem::GridFunction>(&scalar_fespace);
    if (has_d)
        d_gf = std::make_unique<mfem::GridFunction>(&vector_fespace);

    std::filesystem::create_directories(std::filesystem::path(txt_path).parent_path());
    std::ofstream out(txt_path);
    if (!out)
        throw std::runtime_error("Failed to open output txt file: " + txt_path);
    out << std::setprecision(17);

    // Header
    out << "% Model:              FEM_mfem\n"
        << "% Dimension:          " << dim << "\n"
        << "% Nodes:              " << nv << "\n"
        << "% Time steps:         " << nsnaps << "\n"
        << "% Description:        Electric potential, Temperature, Displacement magnitude\n"
        << "% x                       y                        z                        ";
    for (int s = 0; s < nsnaps; ++s)
    {
        if (has_v)
            out << "V @ t=" << times[s] << "  ";
        if (has_t)
            out << "T @ t=" << times[s] << "  ";
        if (has_d)
            out << "disp @ t=" << times[s] << "  ";
    }
    out << "\n";

    // Data — iterate snapshot-major to minimize GF reassignment
    mfem::Vector vec(dim);
    mfem::Vector pt(dim);

    // Pre-compute vertex coordinates
    struct VertexCoord
    {
        double x, y, z;
    };
    std::vector<VertexCoord> coords(nv);
    for (int v = 0; v < nv; ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
            continue;
        auto *tr = mesh.GetElementTransformation(loc.element_id);
        auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
        tr->SetIntPoint(&ip);
        tr->Transform(ip, pt);
        coords[v] = {pt(0), pt(1), dim > 2 ? pt(2) : 0.0};
    }

    // Pre-evaluate all field values: [snap][field_values_per_vertex]
    // Store as flat vectors to avoid per-vertex GF reassignment
    const int nfields = fields_per_snap;
    std::vector<std::vector<double>> snap_data(nsnaps);
    for (int s = 0; s < nsnaps; ++s)
    {
        snap_data[s].resize(nv * nfields);
        if (has_v)
            *v_gf = voltages[s];
        if (has_t)
            *t_gf = temperatures[s];
        if (has_d)
            *d_gf = displacements[s];

        for (int v = 0; v < nv; ++v)
        {
            const auto &loc = locations[v];
            if (loc.element_id < 0)
                continue;
            auto *tr = mesh.GetElementTransformation(loc.element_id);
            auto ip = VertexIntegrationPoint(mesh, loc.element_id, loc.local_vertex_id);
            tr->SetIntPoint(&ip);

            int fi = 0;
            if (has_v)
                snap_data[s][v * nfields + fi++] = v_gf->GetValue(*tr, ip);
            if (has_t)
                snap_data[s][v * nfields + fi++] = t_gf->GetValue(*tr, ip);
            if (has_d)
            {
                d_gf->GetVectorValue(*tr, ip, vec);
                double mag2 = 0.0;
                for (int d = 0; d < dim; ++d)
                    mag2 += vec(d) * vec(d);
                snap_data[s][v * nfields + fi++] = std::sqrt(mag2);
            }
        }
    }

    // Write rows
    for (int v = 0; v < nv; ++v)
    {
        const auto &loc = locations[v];
        if (loc.element_id < 0)
            continue;
        out << coords[v].x << " " << coords[v].y << " " << coords[v].z;
        for (int s = 0; s < nsnaps; ++s)
        {
            for (int fi = 0; fi < nfields; ++fi)
                out << " " << snap_data[s][v * nfields + fi];
        }
        out << "\n";
    }
}

}  // namespace fem::output

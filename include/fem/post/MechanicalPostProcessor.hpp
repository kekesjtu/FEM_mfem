#pragma once

#include "mfem.hpp"

#include <string>

namespace fem::post
{
struct LameParameters
{
    mfem::Coefficient &lambda;
    mfem::Coefficient &mu;
};

struct ThermalStrainParameters
{
    mfem::Coefficient &alpha_sec;
    mfem::Coefficient &temperature;
    double reference_temperature = 293.15;
};

class MechanicalPostProcessor
{
  public:
    static void FillVonMisesNodalField(mfem::Mesh &mesh, mfem::GridFunction &displacement,
                                       const LameParameters &lame,
                                       mfem::FiniteElementSpace &scalar_space,
                                       mfem::GridFunction &von_mises,
                                       const ThermalStrainParameters *thermal = nullptr);

    static void ExportElementStressCsv(const std::string &csv_path, mfem::Mesh &mesh,
                                       mfem::GridFunction &displacement, const LameParameters &lame,
                                       const ThermalStrainParameters *thermal = nullptr);
};
}  // namespace fem::post

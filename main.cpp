#include "fem/app/Application.hpp"

#include "mfem.hpp"
#include "spdlog/spdlog.h"

#include <exception>
#include <string>

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();
#endif

    const std::string config_path =
        (argc > 1) ? argv[1] : "configs/busbar_electro_thermal_iteration.json";

    try
    {
        fem::app::Application app;
        return app.Run(config_path);
    }
    catch (const std::exception &ex)
    {
        spdlog::critical("Fatal error: {}", ex.what());
        return 1;
    }
}

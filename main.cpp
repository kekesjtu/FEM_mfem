#include "app/Application.hpp"

#include "mfem.hpp"
#include "spdlog/spdlog.h"

#include <exception>
#include <string>

int main(int argc, char *argv[])
{
    // MPI_THREAD_SINGLE（级别 0）
    // 进程只有一个线程，且该线程同时负责 MPI 通信和计算。
    // MPI_THREAD_FUNNELED（级别 1）
    // 进程可以有多个线程，但只有主线程（即调用 MPI_Init_thread 的那个线程）可以调用 MPI 函数。
    int provided = MPI_THREAD_SINGLE;

    // MPI::Init 返回provided的值，表示MPI实际提供的线程支持级别
    mfem::Mpi::Init(argc, argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED)
    {
        spdlog::critical("MPI implementation does not provide MPI_THREAD_FUNNELED");
        return 1;
    }
    mfem::Hypre::Init();

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

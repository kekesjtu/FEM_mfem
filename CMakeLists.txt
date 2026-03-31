cmake_minimum_required(VERSION 3.16)

# 项目设置
project(FEM_mfem 
    VERSION 1.0.0
    DESCRIPTION "Finite Element Method Solver based on MFEM"
    LANGUAGES CXX C Fortran
)

# 生成 compile_commands.json 用于 IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 编译选项
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -fdiagnostics-color=auto)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-g)
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
endif()

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# 查找OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "OpenMP not found. Parallel assembly will be disabled.")
endif()

# 依赖: spdlog
find_package(spdlog QUIET)
if(NOT spdlog_FOUND)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.15.3
    )
    FetchContent_MakeAvailable(spdlog)
endif()

# 依赖: nlohmann_json
find_package(nlohmann_json QUIET)
if(NOT nlohmann_json_FOUND)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
    )
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# 依赖: muparser
find_package(muparser QUIET)
if(NOT muparser_FOUND)
    set(ENABLE_SAMPLES OFF CACHE BOOL "Disable muparser samples" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "Disable muparser tests" FORCE)
    FetchContent_Declare(
        muparser
        GIT_REPOSITORY https://github.com/beltoforion/muparser.git
        GIT_TAG v2.3.5
    )
    FetchContent_MakeAvailable(muparser)
endif()

# ==============================================================================
# 配置 MFEM 依赖
# 默认使用仓库内的 mfem/ 源码子目录，避免额外安装步骤。
# 如需切换到系统安装版 MFEM，可开启 FEM_USE_SYSTEM_MFEM。
# ==============================================================================
option(FEM_USE_SYSTEM_MFEM "Use externally installed MFEM package" OFF)

set(MFEM_USE_MPI ON CACHE BOOL "Enable MPI support in MFEM" FORCE)
if (MFEM_USE_MPI)
    message(STATUS "Configuring MFEM with MPI support")
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
    set(MFEM_USE_METIS_5 ON)
    set(MFEM_USE_HYPRE ON)
    set(MFEM_USE_SUITESPARSE ON)
    set(HYPRE_INCLUDE_DIRS /usr/include/hypre)
    set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE.so)

    # MUMPS parallel direct solver
    set(MFEM_USE_MUMPS ON CACHE BOOL "Enable MUMPS direct solver" FORCE)
    set(MUMPS_INCLUDE_DIRS "/usr/include" CACHE STRING "MUMPS include path" FORCE)
    set(MUMPS_LIBRARIES
        "/usr/lib/x86_64-linux-gnu/libdmumps.so"
        "/usr/lib/x86_64-linux-gnu/libmumps_common.so"
        "/usr/lib/x86_64-linux-gnu/libpord.so"
        CACHE STRING "MUMPS libraries" FORCE)
    set(MUMPS_REQUIRED_PACKAGES "MPI" "MPI_Fortran" "METIS"
        "ScaLAPACK" "LAPACK" "BLAS"
        CACHE STRING "MUMPS required packages" FORCE)
    set(ScaLAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so"
        CACHE STRING "ScaLAPACK library" FORCE)

    # MKL PARDISO direct solver (serial, operates on SparseMatrix)
    set(MFEM_USE_MKL_PARDISO ON CACHE BOOL "Enable MKL PARDISO direct solver" FORCE)
    set(MKL_PARDISO_INCLUDE_DIRS "/usr/include/mkl" CACHE STRING "MKL include path" FORCE)
    set(MKL_PARDISO_LIBRARIES
        "/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so"
        "/usr/lib/x86_64-linux-gnu/libmkl_sequential.so"
        "/usr/lib/x86_64-linux-gnu/libmkl_core.so"
        CACHE STRING "MKL PARDISO libraries" FORCE)
    set(MKL_LIBRARY_DIR "lib/x86_64-linux-gnu" CACHE STRING "MKL library subdir" FORCE)

    # MKL Cluster PARDISO (parallel, operates on HypreParMatrix)
    set(MFEM_USE_MKL_CPARDISO ON CACHE BOOL "Enable MKL Cluster PARDISO" FORCE)
    set(MKL_CPARDISO_INCLUDE_DIRS "/usr/include/mkl" CACHE STRING "MKL include path" FORCE)
    set(MKL_CPARDISO_LIBRARIES
        "/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so"
        "/usr/lib/x86_64-linux-gnu/libmkl_sequential.so"
        "/usr/lib/x86_64-linux-gnu/libmkl_core.so"
        "/usr/lib/x86_64-linux-gnu/libmkl_blacs_openmpi_lp64.so"
        CACHE STRING "MKL CPARDISO libraries" FORCE)
    set(MKL_MPI_WRAPPER_LIB "mkl_blacs_openmpi_lp64" CACHE STRING "MKL MPI wrapper" FORCE)
endif()

if(FEM_USE_SYSTEM_MFEM)
    # 需要在配置时提供 MFEM_DIR，例如：
    # cmake -DFEM_USE_SYSTEM_MFEM=ON -DMFEM_DIR="/path/to/mfem/lib/cmake/mfem" ..
    find_package(mfem REQUIRED CONFIG)
    set(FEM_MFEM_TARGET mfem)
else()
    # 作为依赖库构建时，关闭不必要目标，减少配置和构建开销。
    set(MFEM_ENABLE_TESTING OFF CACHE BOOL "Disable MFEM tests" FORCE)
    set(MFEM_ENABLE_EXAMPLES OFF CACHE BOOL "Disable MFEM examples" FORCE)
    set(MFEM_ENABLE_MINIAPPS OFF CACHE BOOL "Disable MFEM miniapps" FORCE)

    add_subdirectory(mfem)
    set(FEM_MFEM_TARGET mfem)
endif()

# ==============================================================================
# 源文件收集
# ==============================================================================
# 主程序 + 分层实现
file(GLOB_RECURSE FEM_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
)

set(FEM_SOURCES
    "${CMAKE_SOURCE_DIR}/main.cpp"
    ${FEM_CORE_SOURCES}
)

# 收集头文件 (可选)
file(GLOB_RECURSE FEM_HEADERS CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/include/*.h"
    "${CMAKE_SOURCE_DIR}/include/*.hpp"
)

# 创建可执行文件
add_executable(my_fem ${FEM_SOURCES} ${FEM_HEADERS})

target_include_directories(my_fem PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

# 链接 MFEM 库
target_link_libraries(my_fem PRIVATE
    ${FEM_MFEM_TARGET}
    spdlog::spdlog
    nlohmann_json::nlohmann_json
    muparser
)

# 链接OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(my_fem PRIVATE OpenMP::OpenMP_CXX)
endif()

# 设置目标属性，输出名称也参考旧项目设为 FEM_solution
set_target_properties(my_fem PROPERTIES
    OUTPUT_NAME "FEM_solution"
)

# 打印一些配置信息
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message(STATUS "Binary directory: ${CMAKE_BINARY_DIR}")

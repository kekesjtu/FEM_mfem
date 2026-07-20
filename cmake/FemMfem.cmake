include_guard(GLOBAL)

option(FEM_USE_SYSTEM_MFEM "Use an externally installed MFEM package" OFF)
option(FEM_ENABLE_MUMPS "Build the MUMPS linear-solver backend" ON)
option(FEM_ENABLE_CPARDISO "Build the MKL Cluster PARDISO backend" ON)

if(FEM_ENABLE_MUMPS)
    enable_language(Fortran)
    find_package(MPI REQUIRED COMPONENTS C CXX Fortran)
else()
    find_package(MPI REQUIRED COMPONENTS C CXX)
endif()

if(FEM_USE_SYSTEM_MFEM)
    # MFEM installs MFEMConfig.cmake, so the package name is case-sensitive.
    # Set MFEM_DIR to the directory containing that file when it is not in a
    # standard CMake prefix.
    find_package(MFEM CONFIG REQUIRED)

    set(_required_mfem_features
        MFEM_USE_MPI
        MFEM_THREAD_SAFE
        MFEM_USE_OPENMP
    )
    if(FEM_ENABLE_MUMPS)
        list(APPEND _required_mfem_features MFEM_USE_MUMPS)
    endif()
    if(FEM_ENABLE_CPARDISO)
        list(APPEND _required_mfem_features MFEM_USE_MKL_CPARDISO)
    endif()

    foreach(_mfem_feature IN LISTS _required_mfem_features)
        if(NOT DEFINED ${_mfem_feature} OR NOT ${${_mfem_feature}})
            message(FATAL_ERROR
                "Installed MFEM is missing required feature ${_mfem_feature}")
        endif()
    endforeach()

    set(FEM_MFEM_TARGET mfem)
    return()
endif()

message(STATUS
    "Configuring bundled MFEM (MUMPS=${FEM_ENABLE_MUMPS}, CPARDISO=${FEM_ENABLE_CPARDISO})")

# MFEM feature selection. These values must be set before add_subdirectory().
set(MFEM_USE_MPI ON CACHE BOOL "Enable MPI support in MFEM" FORCE)
set(MFEM_USE_METIS ON CACHE BOOL "Enable METIS support in MFEM" FORCE)
set(MFEM_USE_SUITESPARSE OFF CACHE BOOL "Disable unused SuiteSparse support" FORCE)
set(MFEM_USE_MUMPS ${FEM_ENABLE_MUMPS} CACHE BOOL "Enable MUMPS direct solver" FORCE)
set(MFEM_USE_MKL_PARDISO OFF CACHE BOOL "Disable unused serial MKL PARDISO" FORCE)
set(MFEM_USE_MKL_CPARDISO ${FEM_ENABLE_CPARDISO}
    CACHE BOOL "Enable MKL Cluster PARDISO" FORCE)
set(MFEM_USE_OPENMP ON CACHE BOOL "Enable MFEM OpenMP backend" FORCE)
set(MFEM_THREAD_SAFE ON CACHE BOOL "Enable thread-safe MFEM kernels" FORCE)

set(MFEM_ENABLE_TESTING OFF CACHE BOOL "Disable MFEM tests" FORCE)
set(MFEM_ENABLE_EXAMPLES OFF CACHE BOOL "Disable MFEM examples" FORCE)
set(MFEM_ENABLE_MINIAPPS OFF CACHE BOOL "Disable MFEM miniapps" FORCE)

# Ubuntu installs HYPRE headers below /usr/include/hypre. Use CMake discovery
# instead of architecture-specific absolute paths.
find_path(FEM_HYPRE_INCLUDE_DIR
    NAMES HYPRE.h
    PATH_SUFFIXES hypre
    REQUIRED
)
find_library(FEM_HYPRE_LIBRARY NAMES HYPRE REQUIRED)
set(HYPRE_INCLUDE_DIRS "${FEM_HYPRE_INCLUDE_DIR}"
    CACHE STRING "HYPRE include directories" FORCE)
set(HYPRE_LIBRARIES "${FEM_HYPRE_LIBRARY}"
    CACHE STRING "HYPRE libraries" FORCE)

# METIS is installed in standard multiarch locations on Ubuntu. Explicitly set
# the detected paths because MFEM's finder otherwise starts with source-tree
# sibling-directory hints.
find_path(FEM_METIS_INCLUDE_DIR NAMES metis.h REQUIRED)
find_library(FEM_METIS_LIBRARY NAMES metis REQUIRED)
set(METIS_INCLUDE_DIRS "${FEM_METIS_INCLUDE_DIR}"
    CACHE STRING "METIS include directories" FORCE)
set(METIS_LIBRARIES "${FEM_METIS_LIBRARY}"
    CACHE STRING "METIS libraries" FORCE)
set(METIS_VERSION_5 ON CACHE BOOL "Use the METIS 5 API" FORCE)
set(MFEM_USE_METIS_5 ON CACHE BOOL "Use the METIS 5 API" FORCE)

if(FEM_ENABLE_MUMPS)
    # MUMPS and its distributed dense linear-algebra dependency.
    find_path(FEM_MUMPS_INCLUDE_DIR NAMES dmumps_c.h REQUIRED)
    find_library(FEM_DMUMPS_LIBRARY NAMES dmumps REQUIRED)
    find_library(FEM_MUMPS_COMMON_LIBRARY NAMES mumps_common REQUIRED)
    find_library(FEM_PORD_LIBRARY NAMES pord REQUIRED)
    find_library(FEM_SCALAPACK_LIBRARY NAMES scalapack-openmpi scalapack REQUIRED)

    set(MUMPS_INCLUDE_DIRS "${FEM_MUMPS_INCLUDE_DIR}"
        CACHE STRING "MUMPS include directories" FORCE)
    set(MUMPS_LIBRARIES
        "${FEM_DMUMPS_LIBRARY};${FEM_MUMPS_COMMON_LIBRARY};${FEM_PORD_LIBRARY}"
        CACHE STRING "MUMPS libraries" FORCE)
    set(MUMPS_REQUIRED_PACKAGES
        "MPI;MPI_Fortran;METIS;ScaLAPACK;LAPACK;BLAS"
        CACHE STRING "MUMPS required packages" FORCE)
    set(ScaLAPACK_LIBRARIES "${FEM_SCALAPACK_LIBRARY}"
        CACHE STRING "ScaLAPACK libraries" FORCE)
endif()

if(FEM_ENABLE_CPARDISO)
    # Ubuntu's MKL packages place headers in include/mkl and libraries in the
    # platform multiarch directory. Let find_path/find_library resolve both.
    find_path(FEM_MKL_INCLUDE_DIR
        NAMES mkl_cluster_sparse_solver.h
        PATH_SUFFIXES mkl
        REQUIRED
    )
    find_library(FEM_MKL_LP64_LIBRARY NAMES mkl_intel_lp64 REQUIRED)
    find_library(FEM_MKL_SEQUENTIAL_LIBRARY NAMES mkl_sequential REQUIRED)
    find_library(FEM_MKL_CORE_LIBRARY NAMES mkl_core REQUIRED)
    find_library(FEM_MKL_BLACS_LIBRARY NAMES mkl_blacs_openmpi_lp64 REQUIRED)
    get_filename_component(FEM_MKL_LIBRARY_DIR "${FEM_MKL_CORE_LIBRARY}" DIRECTORY)

    set(MKL_LIBRARY_DIR "${FEM_MKL_LIBRARY_DIR}"
        CACHE STRING "MKL library directory" FORCE)
    set(MKL_MPI_WRAPPER_LIB "mkl_blacs_openmpi_lp64"
        CACHE STRING "MKL MPI wrapper library" FORCE)
    set(MKL_CPARDISO_INCLUDE_DIRS "${FEM_MKL_INCLUDE_DIR}"
        CACHE STRING "MKL Cluster PARDISO include directories" FORCE)
    set(MKL_CPARDISO_LIBRARIES
        "${FEM_MKL_LP64_LIBRARY};${FEM_MKL_SEQUENTIAL_LIBRARY};${FEM_MKL_CORE_LIBRARY};${FEM_MKL_BLACS_LIBRARY}"
        CACHE STRING "MKL Cluster PARDISO libraries" FORCE)
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/mfem" "${PROJECT_BINARY_DIR}/mfem")
set(FEM_MFEM_TARGET mfem)

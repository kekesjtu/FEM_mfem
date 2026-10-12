include_guard(GLOBAL)

option(FEM_USE_SYSTEM_MFEM "Use an externally installed MFEM package" OFF)
option(FEM_ENABLE_MUMPS "Build the MUMPS linear-solver backend" ON)
option(FEM_ENABLE_CPARDISO "Build the MKL Cluster PARDISO backend" ON)

find_package(MPI REQUIRED COMPONENTS C CXX)

if(FEM_USE_SYSTEM_MFEM)
    # MFEM exports can refer to modern hypre's imported target.
    find_package(HYPRE CONFIG QUIET)
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

if(FEM_ENABLE_MUMPS)
    include(CheckLanguage)
    check_language(Fortran)
    if(NOT CMAKE_Fortran_COMPILER)
        message(FATAL_ERROR
            "Bundled MFEM with MUMPS needs a Fortran compiler and matching MPI Fortran libraries. "
            "On Windows configure Intel ifx in its developer environment, or use "
            "-DFEM_ENABLE_MUMPS=OFF for the MPI/hypre iterative-solver build.")
    endif()
    enable_language(Fortran)
    find_package(MPI REQUIRED COMPONENTS Fortran)
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

# Accept package-specific prefixes, a shared dependency root, and standard
# system locations. Explicit *_INCLUDE_DIRS / *_LIBRARIES remain overridable.
set(HYPRE_ROOT "$ENV{HYPRE_ROOT}" CACHE PATH "HYPRE install prefix")
set(METIS_ROOT "$ENV{METIS_ROOT}" CACHE PATH "METIS install prefix")
set(MUMPS_ROOT "$ENV{MUMPS_ROOT}" CACHE PATH "MUMPS install prefix")
set(MKL_ROOT "$ENV{MKLROOT}" CACHE PATH "oneMKL install prefix")

# Modern hypre exports a target carrying its MPI/BLAS/OpenMP dependencies.
# Ubuntu's older package can still be discovered as headers + library.
find_package(HYPRE CONFIG QUIET)
if(TARGET HYPRE::HYPRE)
    get_target_property(_hypre_includes HYPRE::HYPRE INTERFACE_INCLUDE_DIRECTORIES)
    set(HYPRE_INCLUDE_DIRS "${_hypre_includes}" CACHE STRING "HYPRE includes")
    set(HYPRE_LIBRARIES HYPRE::HYPRE CACHE STRING "HYPRE libraries")
    # MFEM expects an integer, whereas ConfigVersion packages report X.Y.Z.
    if(HYPRE_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)$")
        math(EXPR _hypre_version
            "10000 * ${CMAKE_MATCH_1} + 100 * ${CMAKE_MATCH_2} + ${CMAKE_MATCH_3}")
        set(HYPRE_VERSION "${_hypre_version}")
    endif()
else()
    # A failed Config lookup sets a normal HYPRE_FOUND=FALSE variable, which
    # would shadow the cache result produced later by MFEM's module finder.
    unset(HYPRE_FOUND)
    if(NOT HYPRE_INCLUDE_DIRS)
        find_path(FEM_HYPRE_INCLUDE_DIR NAMES HYPRE.h
            HINTS "${HYPRE_ROOT}" "${FEM_DEPENDENCY_ROOT}/hypre"
            PATH_SUFFIXES include include/hypre hypre REQUIRED)
        set(HYPRE_INCLUDE_DIRS "${FEM_HYPRE_INCLUDE_DIR}" CACHE STRING "HYPRE includes")
    endif()
    if(NOT HYPRE_LIBRARIES)
        find_library(FEM_HYPRE_LIBRARY NAMES HYPRE
            HINTS "${HYPRE_ROOT}" "${FEM_DEPENDENCY_ROOT}/hypre"
            PATH_SUFFIXES lib lib64 lib/x64 REQUIRED)
        set(HYPRE_LIBRARIES "${FEM_HYPRE_LIBRARY}" CACHE STRING "HYPRE libraries")
    endif()
endif()

# METIS is installed in standard multiarch locations on Ubuntu. Explicitly set
# the detected paths because MFEM's finder otherwise starts with source-tree
# sibling-directory hints.
if(NOT METIS_INCLUDE_DIRS)
    find_path(FEM_METIS_INCLUDE_DIR NAMES metis.h
        HINTS "${METIS_ROOT}" "${FEM_DEPENDENCY_ROOT}/metis"
        PATH_SUFFIXES include REQUIRED)
    set(METIS_INCLUDE_DIRS "${FEM_METIS_INCLUDE_DIR}" CACHE STRING "METIS includes")
endif()
if(NOT METIS_LIBRARIES)
    find_library(FEM_METIS_LIBRARY NAMES metis
        HINTS "${METIS_ROOT}" "${FEM_DEPENDENCY_ROOT}/metis"
        PATH_SUFFIXES lib lib64 lib/x64 REQUIRED)
    set(METIS_LIBRARIES "${FEM_METIS_LIBRARY}" CACHE STRING "METIS libraries")
endif()
set(METIS_VERSION_5 ON CACHE BOOL "Use the METIS 5 API" FORCE)
set(MFEM_USE_METIS_5 ON CACHE BOOL "Use the METIS 5 API" FORCE)

if(FEM_ENABLE_MUMPS)
    if(NOT MUMPS_INCLUDE_DIRS)
        find_path(FEM_MUMPS_INCLUDE_DIR NAMES dmumps_c.h
            HINTS "${MUMPS_ROOT}" "${FEM_DEPENDENCY_ROOT}/mumps"
            PATH_SUFFIXES include REQUIRED)
        set(MUMPS_INCLUDE_DIRS "${FEM_MUMPS_INCLUDE_DIR}" CACHE STRING "MUMPS includes")
    endif()
    if(NOT MUMPS_LIBRARIES)
        foreach(_component IN ITEMS DMUMPS MUMPS_COMMON PORD)
            string(TOLOWER "${_component}" _library)
            find_library(FEM_${_component}_LIBRARY NAMES "${_library}"
                HINTS "${MUMPS_ROOT}" "${FEM_DEPENDENCY_ROOT}/mumps"
                PATH_SUFFIXES lib lib64 lib/x64 REQUIRED)
        endforeach()
        if(NOT ScaLAPACK_LIBRARIES)
            find_library(FEM_SCALAPACK_LIBRARY NAMES scalapack-openmpi scalapack
                HINTS "${MUMPS_ROOT}" "${FEM_DEPENDENCY_ROOT}/mumps"
                PATH_SUFFIXES lib lib64 lib/x64 REQUIRED)
            set(ScaLAPACK_LIBRARIES "${FEM_SCALAPACK_LIBRARY}"
                CACHE STRING "ScaLAPACK libraries (include BLACS when using oneMKL)")
        endif()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        # MFEM accepts prefilled MUMPS_LIBRARIES without searching its required
        # packages. Include transitive dependencies explicitly for static builds.
        set(MUMPS_LIBRARIES
            "${FEM_DMUMPS_LIBRARY};${FEM_MUMPS_COMMON_LIBRARY};${FEM_PORD_LIBRARY};${ScaLAPACK_LIBRARIES};${METIS_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${MPI_Fortran_LIBRARIES}"
            CACHE STRING "Complete MUMPS link dependencies")
    endif()
endif()

if(FEM_ENABLE_CPARDISO)
    if(NOT MKL_CPARDISO_INCLUDE_DIRS)
        find_path(FEM_MKL_INCLUDE_DIR NAMES mkl_cluster_sparse_solver.h
            HINTS "${MKL_ROOT}"
            PATH_SUFFIXES include include/mkl mkl REQUIRED)
        set(MKL_CPARDISO_INCLUDE_DIRS "${FEM_MKL_INCLUDE_DIR}"
            CACHE STRING "MKL Cluster PARDISO include directories")
    endif()
    if(WIN32)
        set(_default_blacs mkl_blacs_msmpi_lp64)
    else()
        set(_default_blacs mkl_blacs_openmpi_lp64)
    endif()
    set(MKL_MPI_WRAPPER_LIB "${_default_blacs}" CACHE STRING
        "oneMKL BLACS library matching the MPI implementation")
    if(NOT MKL_CPARDISO_LIBRARIES)
        foreach(_component IN ITEMS LP64 SEQUENTIAL CORE BLACS)
            if(_component STREQUAL "LP64")
                set(_library mkl_intel_lp64)
            elseif(_component STREQUAL "BLACS")
                set(_library "${MKL_MPI_WRAPPER_LIB}")
            else()
                string(TOLOWER "mkl_${_component}" _library)
            endif()
            find_library(FEM_MKL_${_component}_LIBRARY NAMES "${_library}"
                HINTS "${MKL_ROOT}"
                PATH_SUFFIXES lib lib64 lib/intel64 lib/intel64_win REQUIRED)
        endforeach()
        get_filename_component(FEM_MKL_LIBRARY_DIR "${FEM_MKL_CORE_LIBRARY}" DIRECTORY)

        set(MKL_LIBRARY_DIR "${FEM_MKL_LIBRARY_DIR}"
            CACHE STRING "MKL library directory")
        set(MKL_CPARDISO_LIBRARIES
            "${FEM_MKL_LP64_LIBRARY};${FEM_MKL_SEQUENTIAL_LIBRARY};${FEM_MKL_CORE_LIBRARY};${FEM_MKL_BLACS_LIBRARY}"
            CACHE STRING "MKL Cluster PARDISO libraries")
    endif()
endif()

add_subdirectory("${PROJECT_SOURCE_DIR}/mfem" "${PROJECT_BINARY_DIR}/mfem")
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # MFEM contains UTF-8 comments. Decoding them with a legacy Windows code
    # page can consume line endings and hide subsequent source statements.
    target_compile_options(mfem PRIVATE /utf-8)
endif()
if(FEM_ENABLE_MUMPS)
    # Static MUMPS contains Fortran objects although MFEM itself is C++.
    target_link_directories(mfem PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    target_link_libraries(mfem PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()
set(FEM_MFEM_TARGET mfem)

cmake_minimum_required(VERSION 3.23)

project(FEM_mfem
    VERSION 1.0.0
    DESCRIPTION "Finite Element Method Solver based on MFEM"
    LANGUAGES C CXX Fortran
)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(FemProjectSettings)
include(FemDependencies)
include(FemMfem)
include(FemTarget)

fem_print_configuration()

include_guard(GLOBAL)

file(GLOB_RECURSE FEM_CORE_SOURCES CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
)

file(GLOB_RECURSE FEM_HEADERS CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
)

add_executable(my_fem
    "${PROJECT_SOURCE_DIR}/main.cpp"
    ${FEM_CORE_SOURCES}
    ${FEM_HEADERS}
)

target_include_directories(my_fem PRIVATE
    "${PROJECT_SOURCE_DIR}/include"
)

target_link_libraries(my_fem PRIVATE
    fem_project_options
    ${FEM_MFEM_TARGET}
    MPI::MPI_CXX
    OpenMP::OpenMP_CXX
    spdlog::spdlog
    nlohmann_json::nlohmann_json
    muparser::muparser
)

set_target_properties(my_fem PROPERTIES
    OUTPUT_NAME "FEM_solution"
)

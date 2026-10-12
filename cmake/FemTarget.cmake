include_guard(GLOBAL)

add_library(fem_core STATIC)
add_library(fem::core ALIAS fem_core)

target_sources(fem_core PRIVATE
    "${PROJECT_SOURCE_DIR}/src/BCbuilder/MechanicalBCBuilder.cpp"
    "${PROJECT_SOURCE_DIR}/src/BCbuilder/ScalarBCBuilder.cpp"
    "${PROJECT_SOURCE_DIR}/src/app/Application.cpp"
    "${PROJECT_SOURCE_DIR}/src/app/PythonTool.cpp"
    "${PROJECT_SOURCE_DIR}/src/assembler/ElectrostaticAssembler.cpp"
    "${PROJECT_SOURCE_DIR}/src/assembler/LocalAssemblyExecutor.cpp"
    "${PROJECT_SOURCE_DIR}/src/assembler/MechanicalAssembler.cpp"
    "${PROJECT_SOURCE_DIR}/src/assembler/ThermalAssembler.cpp"
    "${PROJECT_SOURCE_DIR}/src/coeff/CoefficientManager.cpp"
    "${PROJECT_SOURCE_DIR}/src/coeff/ElectrostaticCoeffs.cpp"
    "${PROJECT_SOURCE_DIR}/src/coeff/MechanicalCoeffs.cpp"
    "${PROJECT_SOURCE_DIR}/src/coeff/ThermalCoeffs.cpp"
    "${PROJECT_SOURCE_DIR}/src/coupling/PicardCoupler.cpp"
    "${PROJECT_SOURCE_DIR}/src/frontend/ConfigLoader.cpp"
    "${PROJECT_SOURCE_DIR}/src/frontend/Expression.cpp"
    "${PROJECT_SOURCE_DIR}/src/integrator/MechanicalIntegrators.cpp"
    "${PROJECT_SOURCE_DIR}/src/integrator/PoissonIntegrators.cpp"
    "${PROJECT_SOURCE_DIR}/src/log/Logger.cpp"
    "${PROJECT_SOURCE_DIR}/src/output/ParaviewExporter.cpp"
    "${PROJECT_SOURCE_DIR}/src/output/ResultExporter.cpp"
    "${PROJECT_SOURCE_DIR}/src/output/SolutionTextExporter.cpp"
    "${PROJECT_SOURCE_DIR}/src/physics/ElectrostaticFieldSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/physics/MechanicalFieldSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/physics/ThermalFieldSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/solver/LinearSolverFactory.cpp"
    "${PROJECT_SOURCE_DIR}/src/solver/MfemPcgSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/steady/SteadySolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/transient/TransientSolver.cpp"
    "${PROJECT_SOURCE_DIR}/src/transient/TransientUtils.cpp"
)

if(FEM_ENABLE_MUMPS)
    target_sources(fem_core PRIVATE
        "${PROJECT_SOURCE_DIR}/src/solver/MfemMumpsSolver.cpp"
    )
    target_compile_definitions(fem_core PRIVATE FEM_ENABLE_MUMPS)
endif()

if(FEM_ENABLE_CPARDISO)
    target_sources(fem_core PRIVATE
        "${PROJECT_SOURCE_DIR}/src/solver/MfemPardisoSolver.cpp"
    )
    target_compile_definitions(fem_core PRIVATE FEM_ENABLE_CPARDISO)
endif()

target_include_directories(fem_core PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

target_compile_features(fem_core PUBLIC cxx_std_17)

target_link_libraries(fem_core
    PUBLIC
        ${FEM_MFEM_TARGET}
        MPI::MPI_CXX
        OpenMP::OpenMP_CXX
        spdlog::spdlog
        nlohmann_json::nlohmann_json
        muparser::muparser
    PRIVATE
        fem_project_options
)

add_executable(FEM_solution
    "${PROJECT_SOURCE_DIR}/main.cpp"
)

target_link_libraries(FEM_solution PRIVATE
    fem::core
    fem_project_options
)

# Small standalone regression target; no MFEM/MPI linkage and no default build cost.
add_executable(fem_python_tool_tests EXCLUDE_FROM_ALL
    "${PROJECT_SOURCE_DIR}/tests/PythonToolTests.cpp"
    "${PROJECT_SOURCE_DIR}/src/app/PythonTool.cpp"
)
target_include_directories(fem_python_tool_tests PRIVATE "${PROJECT_SOURCE_DIR}/include")
target_compile_features(fem_python_tool_tests PRIVATE cxx_std_17)
target_link_libraries(fem_python_tool_tests PRIVATE fem_project_options)

cmake_minimum_required(VERSION 3.20)
project(cfd_solvers VERSION 0.9.0 LANGUAGES CXX)

option(CFD_ENABLE_OPENMP "Enable multicore CPU execution through OpenMP" ON)
option(CFD_ENABLE_SYCL "Enable AdaptiveCpp/SYCL accelerator backend" OFF)
option(CFD_ENABLE_MPI "Enable MPI distributed-memory runtime" OFF)
option(CFD_ENABLE_NATIVE_ARCH "Enable -march=native on supported compilers" ON)
option(CFD_FAST_MATH "Enable aggressive floating-point optimizations" OFF)
option(CFD_BUILD_TESTS "Build tests" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(cfd_solvers
    src/em/frequency_domain.cpp
    src/em/edge_fem2d.cpp
    src/particle/electromagnetic.cpp
    src/rf/network.cpp
    src/rf/nport.cpp
    src/rf/peec.cpp
    src/rf/antenna.cpp
    src/rf/thin_wire_mom.cpp
    src/circuit/spice.cpp
    src/circuit/analysis.cpp
    src/circuit/raw.cpp
    src/circuit/osdi_loader.cpp
    src/core/task_graph.cpp
    src/physics/material.cpp
    src/physics/material_bridge.cpp
    src/core/mixed_precision.cpp
    src/core/geometric_multigrid.cpp
    src/core/distributed_csr.cpp
    src/workflow/schema.cpp
    src/workflow/provenance.cpp
    src/workflow/validation_catalog.cpp
    src/workflow/runner.cpp
    src/fvm/derived_fields.cpp
    src/fvm/reordering.cpp
    src/fvm/fv_matrix.cpp
    src/fvm/turbulence.cpp
    src/fvm/vof.cpp
    src/fvm/particles.cpp
    src/fvm/compressible1d.cpp
    src/fvm/advanced_models.cpp
    src/fvm/reacting.cpp
    src/fvm/radiation.cpp
    src/fvm/premixed_flame1d.cpp
    src/fvm/radiation_network.cpp
    src/fvm/mesh_motion.cpp
    src/lbm/diagnostics.cpp
    src/lbm/visualization.cpp
    src/lbm/advanced_d2q9.cpp
    src/lbm/voxelize.cpp
    src/fvm/thermal_transport.cpp
    src/fvm/species_transport.cpp
    src/fem/dynamic_bar1d.cpp
    src/fem/stokes2d.cpp
    src/fem/nonlinear_geometry.cpp
    src/fem/brinkman2d.cpp
    src/optics/diffraction.cpp
    src/optics/nonsequential.cpp
    src/core/backend.cpp
    src/core/profiler.cpp
    src/fvm/runtime_schemes.cpp
    src/multiphysics/projection.cpp
    src/multiphysics/cht.cpp
    src/multiphysics/bioheat.cpp
    src/multiphysics/electrochemical_heat.cpp
    src/io/mesh_io.cpp
    src/io/openfoam_io.cpp
    src/optics/transform.cpp
    src/optics/analysis.cpp
    src/optics/wavefront.cpp
    src/optics/optimization.cpp
    src/optics/freeform.cpp
    src/optics/serialization.cpp
    src/electrochemistry/corrosion_models.cpp
    src/lbm/d2q9.cpp
    src/lbm/esoteric_pull.cpp
    src/fvm/diffusion2d.cpp
    src/fvm/poly_mesh.cpp
    src/fvm/operators.cpp
    src/fvm/schemes.cpp
    src/fvm/pressure_velocity.cpp
    src/fvm/projection2d.cpp
    src/fvm/collocated_incompressible.cpp
    src/fvm/incompressible2d.cpp
    src/fvm/scalar_transport.cpp
    src/fem/poisson1d.cpp
    src/fem/mesh2d.cpp
    src/fem/mesh3d.cpp
    src/fem/reference_element.cpp
    src/fem/assembly.cpp
    src/fem/scalar_diffusion2d.cpp
    src/fem/electromagnetics2d.cpp
    src/fem/poisson2d.cpp
    src/fem/poisson3d.cpp
    src/fem/heat2d.cpp
    src/fem/elasticity2d.cpp
    src/fem/axisymmetric_elasticity.cpp
    src/fem/nonlinear_poisson2d.cpp
    src/fem/adaptivity.cpp
    src/fem/darcy2d.cpp
    src/fem/magnetostatics2d.cpp
    src/fem/eddy_current2d.cpp
    src/fem/material_laws.cpp
    src/fem/phase_change.cpp
    src/fem/modal_bar1d.cpp
    src/fdtd/maxwell1d.cpp
    src/fdtd/monitor.cpp
    src/fdtd/port.cpp
    src/fdtd/maxwell3d.cpp
    src/fdtd/vtk_io.cpp
    src/fdtd/postprocess.cpp
    src/optics/ray.cpp
    src/optics/sequential.cpp
    src/optics/materials.cpp
    src/optics/paraxial.cpp
    src/optics/polarization.cpp
    src/optics/birefringence.cpp
    src/optics/gaussian_beam.cpp
    src/chemistry/kinetics.cpp
    src/chemistry/thermo.cpp
    src/chemistry/advanced_thermo.cpp
    src/chemistry/implicit_reactor.cpp
    src/chemistry/reactor_network.cpp
    src/chemistry/aqueous_equilibrium.cpp
    src/chemistry/equilibrium.cpp
    src/chemistry/reaction_path.cpp
    src/electrochemistry/electrochemistry.cpp
    src/electrochemistry/corrosion1d.cpp
    src/electrochemistry/nernst_planck1d.cpp
    src/electrochemistry/nernst_planck_poly.cpp
    src/electrochemistry/mixed_potential.cpp
    src/multiphysics/field_registry.cpp
    src/multiphysics/transfer.cpp
    src/multiphysics/partitioned.cpp
    src/multiphysics/electro_thermal.cpp
    src/fdtd/cylindrical_tm.cpp
    src/fvm/function_objects.cpp
    src/fvm/combustion.cpp
    src/multiphysics/deformation_optics.cpp
    src/multiphysics/block_system.cpp
)
add_library(cfd::solvers ALIAS cfd_solvers)

target_include_directories(cfd_solvers PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(cfd_solvers PUBLIC cxx_std_20)
if(CMAKE_DL_LIBS)
    target_link_libraries(cfd_solvers PRIVATE ${CMAKE_DL_LIBS})
endif()


if(MSVC)
    target_compile_options(cfd_solvers PRIVATE /W4 /permissive-)
    if(CFD_FAST_MATH)
        target_compile_options(cfd_solvers PRIVATE /fp:fast)
    endif()
else()
    target_compile_options(cfd_solvers PRIVATE -Wall -Wextra -Wpedantic -Wconversion -Wshadow)
    if(CFD_ENABLE_NATIVE_ARCH)
        target_compile_options(cfd_solvers PRIVATE -march=native)
    endif()
    if(CFD_FAST_MATH)
        target_compile_options(cfd_solvers PRIVATE -ffast-math)
    endif()
endif()

if(CFD_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(cfd_solvers PUBLIC OpenMP::OpenMP_CXX)
        target_compile_definitions(cfd_solvers PUBLIC CFD_HAS_OPENMP=1)
    else()
        message(WARNING "OpenMP requested but not found; building serial CPU backend")
    endif()
endif()


if(CFD_ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS CXX)
    target_sources(cfd_solvers PRIVATE src/distributed/mpi_runtime.cpp)
    target_link_libraries(cfd_solvers PUBLIC MPI::MPI_CXX)
    target_compile_definitions(cfd_solvers PUBLIC CFD_HAS_MPI=1)

    add_executable(cfd-distributed apps/cfd-distributed.cpp)
    target_link_libraries(cfd-distributed PRIVATE cfd::solvers)
endif()

if(CFD_ENABLE_SYCL)
    find_package(AdaptiveCpp REQUIRED)
    target_sources(cfd_solvers PRIVATE
        src/sycl/d2q9_sycl.cpp
        src/sycl/esoteric_pull_sycl.cpp
    )
    target_compile_definitions(cfd_solvers PUBLIC CFD_HAS_SYCL=1)
    add_sycl_to_target(TARGET cfd_solvers SOURCES
        src/sycl/d2q9_sycl.cpp
        src/sycl/esoteric_pull_sycl.cpp
    )
    if(TARGET cfd-distributed)
        add_sycl_to_target(TARGET cfd-distributed SOURCES apps/cfd-distributed.cpp)
    endif()
endif()

add_executable(cfd-solve apps/cfd-solve.cpp apps/continuation_cases.cpp)
target_link_libraries(cfd-solve PRIVATE cfd::solvers)
if(CFD_ENABLE_SYCL)
    add_sycl_to_target(TARGET cfd-solve SOURCES apps/cfd-solve.cpp)
endif()

add_executable(cfd-bench apps/cfd-bench.cpp)
target_link_libraries(cfd-bench PRIVATE cfd::solvers)

add_executable(cfd-halo-plan apps/cfd-halo-plan.cpp)
target_link_libraries(cfd-halo-plan PRIVATE cfd::solvers)

add_executable(cfd-validate apps/cfd-validate.cpp)
target_link_libraries(cfd-validate PRIVATE cfd::solvers)
if(CFD_ENABLE_SYCL)
    add_sycl_to_target(TARGET cfd-bench SOURCES apps/cfd-bench.cpp)
endif()

include(CTest)
if(CFD_BUILD_TESTS AND BUILD_TESTING)
    add_executable(cfd-tests tests/test_main.cpp)
    target_link_libraries(cfd-tests PRIVATE cfd::solvers)
    if(CFD_ENABLE_SYCL)
        add_sycl_to_target(TARGET cfd-tests SOURCES tests/test_main.cpp)
    endif()
    add_test(NAME cfd-tests COMMAND cfd-tests)
    add_executable(cfd-completion-tests tests/test_completion.cpp)
    target_link_libraries(cfd-completion-tests PRIVATE cfd::solvers)
    add_executable(cfd-spice-io-tests tests/test_spice_io.cpp)
    target_link_libraries(cfd-spice-io-tests PRIVATE cfd::solvers)
    add_test(NAME cfd-spice-io-tests COMMAND cfd-spice-io-tests)
    add_executable(cfd-next-analysis-tests tests/test_next_analysis.cpp)
    target_link_libraries(cfd-next-analysis-tests PRIVATE cfd::solvers)
    add_test(NAME cfd-next-analysis-tests COMMAND cfd-next-analysis-tests)
    add_executable(cfd-cst-expansion-tests tests/test_cst_expansion.cpp)
    target_link_libraries(cfd-cst-expansion-tests PRIVATE cfd::solvers)
    add_test(NAME cfd-cst-expansion-tests COMMAND cfd-cst-expansion-tests)
    add_library(cfd-osdi-fixture SHARED tests/osdi_fixture.cpp)
    add_executable(cfd-osdi-tests tests/test_osdi.cpp)
    target_link_libraries(cfd-osdi-tests PRIVATE cfd::solvers)
    add_test(NAME cfd-osdi-tests COMMAND cfd-osdi-tests $<TARGET_FILE:cfd-osdi-fixture>)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(cfd-completion-tests PRIVATE -O1)
    endif()
    add_test(NAME cfd-completion-tests COMMAND cfd-completion-tests)
    add_executable(cfd-sanitize-smoke tests/test_sanitize_smoke.cpp)
    target_link_libraries(cfd-sanitize-smoke PRIVATE cfd::solvers)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Regression drivers are not performance targets. Avoid spending release-build
        # time applying -O3 to large test-only translation units.
        target_compile_options(cfd-tests PRIVATE -O1)
        target_compile_options(cfd-sanitize-smoke PRIVATE -O1)
    endif()
    if(CFD_ENABLE_SYCL)
        add_sycl_to_target(TARGET cfd-sanitize-smoke SOURCES tests/test_sanitize_smoke.cpp)
    endif()
    add_test(NAME cfd-sanitize-smoke COMMAND cfd-sanitize-smoke)
    add_test(NAME cfd-bench-smoke COMMAND cfd-bench --lattice d3q19 --nx 8 --ny 8 --nz 8 --warmup 1 --steps 2 --threads 1 --csv)
    add_test(NAME cfd-bench-json-smoke COMMAND cfd-bench --lattice d3q19 --nx 8 --ny 8 --nz 8 --warmup 1 --steps 2 --threads 1 --json)
    add_test(NAME cfd-halo-plan-smoke COMMAND cfd-halo-plan --lattice d3q27 --nx 8 --ny 8 --nz 8)
    add_test(NAME cfd-validation-catalog-smoke COMMAND cfd-validate)
    foreach(case_name lbm-trt fvm-workspace fdtd-pml1d chemistry-reactor chemistry-equilibrium optics-gaussian multiphysics-thermoelastic em-frequency1d em-edge2d particle-pic1d bioheat-sar rf-dipole rf-microstrip rf-multiwire spice-rc spice-diode spice-adaptive spice-pss-pz)
        add_test(NAME cfd-${case_name}-smoke COMMAND cfd-solve ${case_name})
    endforeach()
    add_test(NAME cfd-fvm-operators-smoke COMMAND cfd-solve fvm-operators3d)
    add_test(NAME cfd-fvm-projection-smoke COMMAND cfd-solve fvm-projection2d)
    add_test(NAME cfd-fvm-navier-stokes-smoke COMMAND cfd-solve fvm-taylor-green2d)
    add_test(NAME cfd-fvm-collocated-skew-smoke COMMAND cfd-solve fvm-collocated-skew)
    add_test(NAME cfd-fvm-collocated-channel-smoke COMMAND cfd-solve fvm-collocated-channel)
    add_test(NAME cfd-fvm-scalar-transport-smoke COMMAND cfd-solve fvm-scalar-transport)
    add_test(NAME cfd-fem-poisson3d-smoke COMMAND cfd-solve fem-poisson3d)
    add_test(NAME cfd-fem-nonlinear-smoke COMMAND cfd-solve fem-nonlinear-poisson)
    add_test(NAME cfd-fem-darcy-smoke COMMAND cfd-solve fem-darcy)
    add_test(NAME cfd-fem-magnetostatic-smoke COMMAND cfd-solve fem-magnetostatic)
    add_test(NAME cfd-fem-modal-smoke COMMAND cfd-solve fem-modal-bar)
    add_test(NAME cfd-fdtd-mur1d-smoke COMMAND cfd-solve fdtd-mur1d)
    add_test(NAME cfd-fdtd-dispersive-smoke COMMAND cfd-solve fdtd-dispersive1d)
    add_test(NAME cfd-fdtd-port-vtk-smoke COMMAND cfd-solve fdtd-port-vtk)
    add_test(NAME cfd-fdtd-cpml-tfsf-smoke COMMAND cfd-solve fdtd-cpml-tfsf)
    add_test(NAME cfd-fdtd-lumped-rlc-smoke COMMAND cfd-solve fdtd-lumped-rlc)
    add_test(NAME cfd-multiphysics-coupling-smoke COMMAND cfd-solve multiphysics-coupling)
    add_test(NAME cfd-multiphysics-electrothermal-smoke COMMAND cfd-solve multiphysics-electrothermal)
    add_test(NAME cfd-optics-lens-smoke COMMAND cfd-solve optics-lens)
    add_test(NAME cfd-electrochem-corrosion-smoke COMMAND cfd-solve electrochem-corrosion1d)
    add_test(NAME cfd-electrochem-pnp-smoke COMMAND cfd-solve electrochem-pnp1d)
    add_test(NAME cfd-electrochem-pnp-poly-smoke COMMAND cfd-solve electrochem-pnp-poly)
    add_test(NAME cfd-electrochem-galvanic-smoke COMMAND cfd-solve electrochem-galvanic)

    if(CFD_ENABLE_MPI)
        add_executable(cfd-mpi-tests tests/test_mpi.cpp)
        target_link_libraries(cfd-mpi-tests PRIVATE cfd::solvers)
        if(CFD_ENABLE_SYCL)
            add_sycl_to_target(TARGET cfd-mpi-tests SOURCES tests/test_mpi.cpp)
        endif()
        add_test(NAME cfd-mpi-tests
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfd-mpi-tests> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(cfd-mpi-tests PROPERTIES PROCESSORS 4)
        add_test(NAME cfd-mpi-cli-d3q19
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfd-distributed>
                    --backend cpu --lattice d3q19 --nx 18 --ny 14 --nz 10 --steps 4
                    ${MPIEXEC_POSTFLAGS})
        add_test(NAME cfd-mpi-cli-d3q27
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:cfd-distributed>
                    --backend cpu --lattice d3q27 --nx 18 --ny 14 --nz 10 --steps 4
                    ${MPIEXEC_POSTFLAGS})
        set_tests_properties(cfd-mpi-cli-d3q19 cfd-mpi-cli-d3q27 PROPERTIES PROCESSORS 4)
    endif()
endif()

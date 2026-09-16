cmake_minimum_required(VERSION 3.20)
project(cfd_solvers VERSION 0.4.2 LANGUAGES CXX)

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
    src/fem/poisson1d.cpp
    src/fdtd/maxwell1d.cpp
    src/optics/ray.cpp
)
add_library(cfd::solvers ALIAS cfd_solvers)

target_include_directories(cfd_solvers PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(cfd_solvers PUBLIC cxx_std_20)

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

add_executable(cfd-solve apps/cfd-solve.cpp)
target_link_libraries(cfd-solve PRIVATE cfd::solvers)
if(CFD_ENABLE_SYCL)
    add_sycl_to_target(TARGET cfd-solve SOURCES apps/cfd-solve.cpp)
endif()

add_executable(cfd-bench apps/cfd-bench.cpp)
target_link_libraries(cfd-bench PRIVATE cfd::solvers)

add_executable(cfd-halo-plan apps/cfd-halo-plan.cpp)
target_link_libraries(cfd-halo-plan PRIVATE cfd::solvers)
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
    add_executable(cfd-sanitize-smoke tests/test_sanitize_smoke.cpp)
    target_link_libraries(cfd-sanitize-smoke PRIVATE cfd::solvers)
    if(CFD_ENABLE_SYCL)
        add_sycl_to_target(TARGET cfd-sanitize-smoke SOURCES tests/test_sanitize_smoke.cpp)
    endif()
    add_test(NAME cfd-sanitize-smoke COMMAND cfd-sanitize-smoke)
    add_test(NAME cfd-bench-smoke COMMAND cfd-bench --lattice d3q19 --nx 8 --ny 8 --nz 8 --warmup 1 --steps 2 --threads 1 --csv)
    add_test(NAME cfd-halo-plan-smoke COMMAND cfd-halo-plan --lattice d3q27 --nx 8 --ny 8 --nz 8)
    add_test(NAME cfd-fvm-operators-smoke COMMAND cfd-solve fvm-operators3d)
    add_test(NAME cfd-fvm-projection-smoke COMMAND cfd-solve fvm-projection2d)
    add_test(NAME cfd-fvm-navier-stokes-smoke COMMAND cfd-solve fvm-taylor-green2d)
    add_test(NAME cfd-fvm-collocated-skew-smoke COMMAND cfd-solve fvm-collocated-skew)
    add_test(NAME cfd-fvm-collocated-channel-smoke COMMAND cfd-solve fvm-collocated-channel)

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

# syntax=docker/dockerfile:1
#
# cfd_solvers build + runtime image.
#
# Multi-stage: a builder stage compiles the C++20 solvers (Release, OpenMP,
# optional MPI) and a slim runtime stage ships the binaries plus the
# cfd_solvers source tree for the Python ctypes package and the examples.
#
# The same image powers local runs, `docker compose` deployment and CI.

ARG DEBIAN_VERSION=12

# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------
FROM debian:${DEBIAN_VERSION} AS builder
ARG CMAKE_BUILD_TYPE=Release
ARG CFD_ENABLE_OPENMP=ON
ARG CFD_ENABLE_MPI=OFF
ARG CFD_ENABLE_NATIVE_ARCH=ON
ARG CFD_ENABLE_HDF5=ON
ARG CFD_BUILD_TESTS=ON
ARG CFD_BUILD_EXAMPLES=ON
ARG CFD_BUILD_PYTHON_BINDINGS=ON

ENV DEBIAN_FRONTEND=noninteractive
ENV CC=gcc CXX=g++
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build g++ \
        libopenmpi-dev openmpi-bin \
        libhdf5-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY CMakeLists.txt ./
COPY include/ ./include/
COPY src/ ./src/
COPY apps/ ./apps/
COPY examples/ ./examples/
COPY tests/ ./tests/

RUN cmake -S . -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
        -DCFD_ENABLE_OPENMP=${CFD_ENABLE_OPENMP} \
        -DCFD_ENABLE_MPI=${CFD_ENABLE_MPI} \
        -DCFD_ENABLE_NATIVE_ARCH=${CFD_ENABLE_NATIVE_ARCH} \
        -DCFD_ENABLE_HDF5=${CFD_ENABLE_HDF5} \
        -DCFD_BUILD_TESTS=${CFD_BUILD_TESTS} \
        -DCFD_BUILD_EXAMPLES=${CFD_BUILD_EXAMPLES} \
        -DCFD_BUILD_PYTHON_BINDINGS=${CFD_BUILD_PYTHON_BINDINGS} \
    && cmake --build build --parallel "$(nproc)"

# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------
FROM debian:${DEBIAN_VERSION} AS runtime

ARG CFD_ENABLE_MPI=OFF
ENV DEBIAN_FRONTEND=noninteractive
# libhdf5-103 is the bookworm runtime library; openmpi-bin only when MPI is on.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-103 \
        python3 \
    && if [ "${CFD_ENABLE_MPI}" = "ON" ]; then \
         apt-get install -y --no-install-recommends openmpi-bin; \
       fi \
    && rm -rf /var/lib/apt/lists/*

# Ship the source tree (Python package + examples) alongside the binaries so the
# image is self-contained for both CLI and Python usage.
COPY --from=builder /src /cfd_solvers
COPY --from=builder /src/build /cfd_solvers/build

# Non-root execution by default.
RUN useradd --create-home --shell /usr/sbin/nologin cfd \
    && chown -R cfd:cfd /cfd_solvers

USER cfd
WORKDIR /cfd_solvers

ENV CFD_BUILD_DIR=/cfd_solvers/build
ENV OMP_PROC_BIND=spread
ENV OMP_DYNAMIC=false

# Convenience entrypoint: `cfd-solve <case>` by default, override for mpirun/ctest.
ENTRYPOINT ["/cfd_solvers/build/cfd-solve"]
CMD ["--list"]

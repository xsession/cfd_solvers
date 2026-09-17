#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/mpi_runtime.hpp"
#include "cfd/distributed/mpi_sparse.hpp"
#include "cfd/multibody/distributed_dem.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/distributed/mpi_selective_halo.hpp"
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include "cfd/solvers/lbm/precision_pull.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/distributed/mpi_sycl_halo.hpp"
#include "cfd/solvers/lbm/distributed_sycl_pull.hpp"
#endif

#if !defined(CFD_HAS_MPI)
#error "test_mpi requires CFD_HAS_MPI"
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

constexpr float tau = 0.73F;
constexpr float amplitude = 0.017F;
constexpr std::size_t steps = 7;

template<class Descriptor>
void gather_and_verify(const cfd::distributed::MpiCartesianRuntime& runtime,
                       const cfd::lbm::MacroscopicFields& local,
                       cfd::distributed::Extent3 global,
                       const cfd::distributed::Brick& b,
                       const char* label) {
    const int local_cells = static_cast<int>(b.extent.cells());
    std::vector<unsigned long long> local_indices(static_cast<std::size_t>(local_cells));
    std::vector<float> local_values(static_cast<std::size_t>(local_cells) * 4U);
    std::size_t n = 0;
    for (std::size_t z = 0; z < b.extent.z; ++z) {
        for (std::size_t y = 0; y < b.extent.y; ++y) {
            for (std::size_t x = 0; x < b.extent.x; ++x, ++n) {
                const std::size_t gx = b.begin.x + x;
                const std::size_t gy = b.begin.y + y;
                const std::size_t gz = b.begin.z + z;
                local_indices[n] = static_cast<unsigned long long>((gz * global.y + gy) * global.x + gx);
                local_values[4U * n + 0U] = local.rho[n];
                local_values[4U * n + 1U] = local.ux[n];
                local_values[4U * n + 2U] = local.uy[n];
                local_values[4U * n + 3U] = local.uz[n];
            }
        }
    }

    std::vector<int> counts;
    if (runtime.rank() == 0) counts.resize(static_cast<std::size_t>(runtime.size()));
    MPI_Gather(&local_cells, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, runtime.communicator());

    std::vector<int> displacements;
    std::vector<int> value_counts;
    std::vector<int> value_displacements;
    std::vector<unsigned long long> all_indices;
    std::vector<float> all_values;
    if (runtime.rank() == 0) {
        displacements.resize(counts.size());
        value_counts.resize(counts.size());
        value_displacements.resize(counts.size());
        int total = 0;
        int value_total = 0;
        for (std::size_t i = 0; i < counts.size(); ++i) {
            displacements[i] = total;
            value_displacements[i] = value_total;
            value_counts[i] = 4 * counts[i];
            total += counts[i];
            value_total += value_counts[i];
        }
        all_indices.resize(static_cast<std::size_t>(total));
        all_values.resize(static_cast<std::size_t>(value_total));
    }

    MPI_Gatherv(local_indices.data(), local_cells, MPI_UNSIGNED_LONG_LONG,
                all_indices.data(), counts.data(), displacements.data(), MPI_UNSIGNED_LONG_LONG,
                0, runtime.communicator());
    const int local_value_count = 4 * local_cells;
    MPI_Gatherv(local_values.data(), local_value_count, MPI_FLOAT,
                all_values.data(), value_counts.data(), value_displacements.data(), MPI_FLOAT,
                0, runtime.communicator());

    int ok = 1;
    if (runtime.rank() == 0) {
        cfd::lbm::PrecisionPullSolver<Descriptor, float> reference({global.x, global.y, global.z, tau});
        reference.initialize_taylor_green(amplitude);
        reference.step(steps);
        const auto serial = reference.compute_macroscopic();
        double max_error = 0.0;
        for (std::size_t i = 0; i < all_indices.size(); ++i) {
            const auto gi = static_cast<std::size_t>(all_indices[i]);
            max_error = std::max(max_error, std::abs(static_cast<double>(all_values[4U * i + 0U] - serial.rho[gi])));
            max_error = std::max(max_error, std::abs(static_cast<double>(all_values[4U * i + 1U] - serial.ux[gi])));
            max_error = std::max(max_error, std::abs(static_cast<double>(all_values[4U * i + 2U] - serial.uy[gi])));
            max_error = std::max(max_error, std::abs(static_cast<double>(all_values[4U * i + 3U] - serial.uz[gi])));
        }
        if (max_error >= 2.0e-6 || all_indices.size() != global.cells()) {
            std::cerr << label << " parity failed: max_error=" << max_error << '\n';
            ok = 0;
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, runtime.communicator());
    if (!ok) throw std::runtime_error(std::string(label) + " parity regression");
}

template<class Descriptor>
void run_cpu_case(const cfd::distributed::MpiCartesianRuntime& runtime) {
    const auto global = runtime.brick().global;
    cfd::lbm::DistributedPullBlock<Descriptor> block({global, tau}, runtime.brick());
    block.initialize_taylor_green(amplitude);
    cfd::distributed::MpiSelectiveHaloExchange<Descriptor> exchange(runtime, block.current_grid());

    for (std::size_t i = 0; i < steps; ++i) {
        exchange.begin();
        block.compute_strict_interior();
        exchange.wait();
        block.compute_boundary();
        block.finish_step();
    }
    gather_and_verify<Descriptor>(runtime, block.compute_macroscopic(), global, block.brick(), "MPI CPU selective");
}

#if defined(CFD_HAS_SYCL)
template<class Descriptor>
void run_sycl_case(const cfd::distributed::MpiCartesianRuntime& runtime) {
    const auto global = runtime.brick().global;
    const auto device = cfd::distributed::device_for_local_rank(
        static_cast<std::size_t>(runtime.local_rank()),
        static_cast<std::size_t>(runtime.local_size()));
    cfd::lbm::DistributedSyclPullBlock<Descriptor> block({global, tau}, runtime.brick(), device);
    block.initialize_taylor_green(amplitude);
    cfd::distributed::MpiSyclSelectiveHaloExchange<Descriptor> exchange(
        runtime, block, cfd::distributed::GpuMpiMode::staged_host);

    for (std::size_t i = 0; i < steps; ++i) {
        exchange.begin();
        block.compute_strict_interior();
        exchange.wait();
        block.compute_boundary();
        block.finish_step();
    }
    block.wait();
    gather_and_verify<Descriptor>(runtime, block.download_macroscopic(), global, block.brick(), "MPI SYCL staged");
}
#endif

void run_distributed_sparse_krylov(const cfd::distributed::MpiCartesianRuntime& runtime) {
    const std::size_t ranks = static_cast<std::size_t>(runtime.size());
    const std::size_t rank = static_cast<std::size_t>(runtime.rank());
    const std::size_t n = std::max<std::size_t>(8U, ranks * 4U);
    const std::size_t begin = n * rank / ranks;
    const std::size_t end = n * (rank + 1U) / ranks;

    cfd::core::CsrBuilder builder(n, n);
    for (std::size_t i = 0U; i < n; ++i) {
        if (i > 0U) builder.add(i, i - 1U, -1.0);
        builder.add(i, i, 2.0);
        if (i + 1U < n) builder.add(i, i + 1U, -1.0);
    }
    const auto matrix = builder.build();
    std::vector<double> exact(n);
    for (std::size_t i = 0U; i < n; ++i) exact[i] = 1.0 + 0.125 * static_cast<double>(i);
    std::vector<double> global_rhs(n);
    matrix.multiply(exact, global_rhs);

    std::vector<double> local_rhs(global_rhs.begin() + static_cast<std::ptrdiff_t>(begin),
                                  global_rhs.begin() + static_cast<std::ptrdiff_t>(end));
    std::vector<double> local_x(end - begin, 0.0);
    cfd::distributed::MpiDistributedCsrOperator distributed(matrix, begin, end, runtime.communicator());
    cfd::core::KrylovWorkspace workspace;
    const auto result = cfd::distributed::mpi_distributed_conjugate_gradient(
        distributed, local_rhs, local_x, workspace, 256U, 1.0e-12);
    if (!result.converged) throw std::runtime_error("MPI distributed CG failed to converge");

    double local_error = 0.0;
    for (std::size_t i = begin; i < end; ++i) {
        local_error = std::max(local_error, std::abs(local_x[i - begin] - exact[i]));
    }
    double global_error = 0.0;
    MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, runtime.communicator());
    if (global_error >= 1.0e-9) throw std::runtime_error("MPI distributed CG solution regression");

    if (runtime.size() > 1 && distributed.halo_value_count() == 0U) {
        throw std::runtime_error("MPI distributed CSR partition unexpectedly has no halo columns");
    }
}


void run_distributed_dem_exchange(const cfd::distributed::MpiCartesianRuntime& runtime) {
    using namespace cfd::multibody;
    DemSlabDecomposition decomposition{0.0,static_cast<double>(runtime.size()),runtime.size()};
    MpiDemDomainExchange exchange(decomposition,runtime.communicator());
    const int destination=(runtime.rank()+1)%runtime.size();
    DistributedDemParticle particle;
    particle.global_id=static_cast<std::uint64_t>(runtime.rank()+1000);
    particle.state.position={static_cast<double>(destination)+0.25,0.0,0.0};
    particle.radius=0.05;
    particle.owner_rank=runtime.rank();
    std::vector<DistributedDemParticle> owned{particle};
    auto ghosts=exchange.exchange(owned,0.30);
    if(owned.size()!=1U || owned.front().owner_rank!=runtime.rank()) throw std::runtime_error("MPI DEM ownership migration regression");
    if(decomposition.owner_rank(owned.front().state.position.x)!=runtime.rank()) throw std::runtime_error("MPI DEM migrated particle on wrong rank");
    const int local_ghosts=static_cast<int>(ghosts.size());
    int global_ghosts=0;
    MPI_Allreduce(&local_ghosts,&global_ghosts,1,MPI_INT,MPI_SUM,runtime.communicator());
    const int expected=std::max(0,runtime.size()-1);
    if(global_ghosts!=expected) throw std::runtime_error("MPI DEM ghost exchange count regression");
}

} // namespace

int main(int argc, char** argv) {
    try {
        cfd::distributed::MpiEnvironment environment(argc, argv);
        cfd::distributed::MpiCartesianRuntime runtime({18, 14, 10}, 3, {true, true, true});
        run_cpu_case<cfd::lbm::D3Q19Descriptor>(runtime);
        run_cpu_case<cfd::lbm::D3Q27Descriptor>(runtime);
        run_distributed_sparse_krylov(runtime);
        run_distributed_dem_exchange(runtime);
#if defined(CFD_HAS_SYCL)
        run_sycl_case<cfd::lbm::D3Q19Descriptor>(runtime);
        run_sycl_case<cfd::lbm::D3Q27Descriptor>(runtime);
#endif
        if (runtime.rank() == 0) std::cout << "MPI distributed tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "MPI test failure: " << e.what() << '\n';
        return 1;
    }
}

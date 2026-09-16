#include "cfd/distributed/device_assignment.hpp"
#include "cfd/distributed/mpi_runtime.hpp"
#include "cfd/distributed/mpi_selective_halo.hpp"
#include "cfd/solvers/lbm/distributed_pull.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/distributed/mpi_sycl_halo.hpp"
#include "cfd/solvers/lbm/distributed_sycl_pull.hpp"
#include <sycl/sycl.hpp>
#endif

#if !defined(CFD_HAS_MPI)
#error "cfd-distributed requires CFD_HAS_MPI"
#endif

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

struct Options {
    std::size_t nx{128};
    std::size_t ny{128};
    std::size_t nz{128};
    std::size_t steps{100};
    float tau{0.7F};
    std::string lattice{"d3q19"};
    std::string backend{"cpu"};
    bool gpu_aware_mpi{false};
    bool csv{false};
    std::string checkpoint_prefix;
    std::string restart_prefix;
};

std::size_t parse_size(std::string_view value, std::string_view name) {
    const unsigned long long parsed = std::stoull(std::string(value));
    if (parsed == 0) throw std::invalid_argument(std::string(name) + " must be positive");
    return static_cast<std::size_t>(parsed);
}

Options parse_options(int argc, char** argv) {
    Options result;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        auto require_value = [&](std::string_view name) -> std::string_view {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(name) + " requires a value");
            return argv[++i];
        };
        if (arg == "--nx") result.nx = parse_size(require_value(arg), arg);
        else if (arg == "--ny") result.ny = parse_size(require_value(arg), arg);
        else if (arg == "--nz") result.nz = parse_size(require_value(arg), arg);
        else if (arg == "--steps") result.steps = parse_size(require_value(arg), arg);
        else if (arg == "--tau") result.tau = std::stof(std::string(require_value(arg)));
        else if (arg == "--lattice") result.lattice = std::string(require_value(arg));
        else if (arg == "--backend") result.backend = std::string(require_value(arg));
        else if (arg == "--gpu-aware-mpi") result.gpu_aware_mpi = true;
        else if (arg == "--csv") result.csv = true;
        else if (arg == "--checkpoint") result.checkpoint_prefix = std::string(require_value(arg));
        else if (arg == "--restart") result.restart_prefix = std::string(require_value(arg));
        else if (arg == "--help") {
            std::cout << "Usage: cfd-distributed [--nx N --ny N --nz N --steps N --tau T] "
                         "[--lattice d3q19|d3q27] [--backend cpu|sycl] [--gpu-aware-mpi] [--csv] "
                         "[--checkpoint PREFIX] [--restart PREFIX]\n\n"
                         "The SYCL backend defaults to pinned-host MPI staging. --gpu-aware-mpi passes "
                         "device USM buffers directly to MPI and must only be used with an accelerator-aware MPI stack.\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + std::string(arg));
        }
    }
    if (result.lattice != "d3q19" && result.lattice != "d3q27") {
        throw std::invalid_argument("--lattice must be d3q19 or d3q27");
    }
    if (result.backend != "cpu" && result.backend != "sycl") {
        throw std::invalid_argument("--backend must be cpu or sycl");
    }
    return result;
}

std::string rank_path(const std::string& prefix, int rank) {
    return prefix + ".rank" + std::to_string(rank) + ".bin";
}

void reduce_and_report(const Options& options,
                       const cfd::distributed::MpiCartesianRuntime& runtime,
                       const cfd::distributed::Extent3 global,
                       std::string_view backend,
                       std::string_view lattice,
                       double local_mass,
                       double seconds,
                       std::size_t local_halo_messages,
                       std::size_t local_halo_bytes) {
    double global_mass = 0.0;
    MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM, runtime.communicator());
    double max_seconds = 0.0;
    MPI_Reduce(&seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, runtime.communicator());
    const unsigned long long local_messages = static_cast<unsigned long long>(local_halo_messages);
    const unsigned long long local_bytes = static_cast<unsigned long long>(local_halo_bytes);
    unsigned long long total_messages = 0;
    unsigned long long total_bytes = 0;
    MPI_Reduce(&local_messages, &total_messages, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, runtime.communicator());
    MPI_Reduce(&local_bytes, &total_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, runtime.communicator());

    if (runtime.rank() == 0) {
        const auto pg = runtime.process_grid();
        const double mlups = static_cast<double>(global.cells()) * static_cast<double>(options.steps) /
                             std::max(max_seconds, 1.0e-30) / 1.0e6;
        const double gib = static_cast<double>(total_bytes) * static_cast<double>(options.steps) /
                           (1024.0 * 1024.0 * 1024.0);
        if (options.csv) {
            std::cout << backend << ',' << lattice << ',' << runtime.size() << ','
                      << pg.x << ',' << pg.y << ',' << pg.z << ','
                      << global.x << ',' << global.y << ',' << global.z << ','
                      << options.steps << ',' << std::fixed << std::setprecision(6) << max_seconds << ','
                      << std::setprecision(3) << mlups << ','
                      << total_messages << ',' << total_bytes << ','
                      << std::setprecision(6) << gib << ','
                      << std::setprecision(9) << global_mass << '\n';
        } else {
            std::cout << "backend=" << backend << " lattice=" << lattice
                      << " ranks=" << runtime.size()
                      << " process_grid=" << pg.x << 'x' << pg.y << 'x' << pg.z
                      << " global=" << global.x << 'x' << global.y << 'x' << global.z
                      << " steps=" << options.steps
                      << " seconds=" << std::fixed << std::setprecision(6) << max_seconds
                      << " mlups=" << std::setprecision(3) << mlups
                      << " halo_messages_per_step=" << total_messages
                      << " halo_bytes_per_step=" << total_bytes
                      << " halo_gib=" << std::setprecision(6) << gib
                      << " mass=" << std::setprecision(9) << global_mass << '\n';
        }
    }
}

template<class Descriptor>
int run_cpu(const Options& options,
            const cfd::distributed::MpiCartesianRuntime& runtime,
            cfd::distributed::Extent3 global,
            std::string_view lattice) {
    const cfd::lbm::DistributedLbmConfig config{global, options.tau};
    cfd::lbm::DistributedPullBlock<Descriptor> block(config, runtime.brick());
    if (!options.restart_prefix.empty()) block.load_checkpoint(rank_path(options.restart_prefix, runtime.rank()));
    else block.initialize_taylor_green(0.02F);

    cfd::distributed::MpiSelectiveHaloExchange<Descriptor> exchange(runtime, block.current_grid());

    MPI_Barrier(runtime.communicator());
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t step = 0; step < options.steps; ++step) {
        exchange.begin();
        block.compute_strict_interior();
        exchange.wait();
        block.compute_boundary();
        block.finish_step();
    }
    MPI_Barrier(runtime.communicator());
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (!options.checkpoint_prefix.empty()) {
        block.save_checkpoint(rank_path(options.checkpoint_prefix, runtime.rank()));
        MPI_Barrier(runtime.communicator());
    }

    reduce_and_report(options, runtime, global, "mpi+openmp-selective", lattice,
                      block.local_mass(), seconds, exchange.message_count(), exchange.bytes_per_step());
    return 0;
}

#if defined(CFD_HAS_SYCL)
template<class Descriptor>
int run_sycl(const Options& options,
             const cfd::distributed::MpiCartesianRuntime& runtime,
             cfd::distributed::Extent3 global,
             std::string_view lattice) {
    const auto accelerators = cfd::distributed::visible_accelerators();
    const auto device = cfd::distributed::device_for_local_rank(
        static_cast<std::size_t>(runtime.local_rank()),
        static_cast<std::size_t>(runtime.local_size()));
    if (!options.csv && !accelerators.empty()) {
        const auto assignment = cfd::distributed::assign_device(
            static_cast<std::size_t>(runtime.local_rank()),
            static_cast<std::size_t>(runtime.local_size()), accelerators.size());
        std::cout << "rank=" << runtime.rank()
                  << " device_ordinal=" << assignment.device_ordinal
                  << " device=\"" << device.get_info<sycl::info::device::name>() << "\""
                  << (assignment.oversubscribed ? " oversubscribed=1" : " oversubscribed=0") << '\n';
    } else if (!options.csv) {
        std::cout << "rank=" << runtime.rank()
                  << " device=\"" << device.get_info<sycl::info::device::name>()
                  << "\" accelerator_fallback=1\n";
    }

    const cfd::lbm::DistributedLbmConfig config{global, options.tau};
    cfd::lbm::DistributedSyclPullBlock<Descriptor> block(config, runtime.brick(), device);
    if (!options.restart_prefix.empty()) block.load_checkpoint(rank_path(options.restart_prefix, runtime.rank()));
    else block.initialize_taylor_green(0.02F);

    const auto mode = options.gpu_aware_mpi ? cfd::distributed::GpuMpiMode::direct_device
                                            : cfd::distributed::GpuMpiMode::staged_host;
    cfd::distributed::MpiSyclSelectiveHaloExchange<Descriptor> exchange(runtime, block, mode);

    MPI_Barrier(runtime.communicator());
    const auto start = std::chrono::steady_clock::now();
    for (std::size_t step = 0; step < options.steps; ++step) {
        exchange.begin();
        block.compute_strict_interior();
        exchange.wait();
        block.compute_boundary();
        block.finish_step();
    }
    block.wait();
    MPI_Barrier(runtime.communicator());
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

    if (!options.checkpoint_prefix.empty()) {
        block.save_checkpoint(rank_path(options.checkpoint_prefix, runtime.rank()));
        MPI_Barrier(runtime.communicator());
    }

    reduce_and_report(options, runtime, global,
                      options.gpu_aware_mpi ? "mpi+sycl-direct" : "mpi+sycl-staged",
                      lattice, block.local_mass(), seconds, exchange.message_count(), exchange.bytes_per_step());
    return 0;
}
#endif

} // namespace

int main(int argc, char** argv) {
    try {
        cfd::distributed::MpiEnvironment environment(argc, argv);
        const auto options = parse_options(argc, argv);
        const cfd::distributed::Extent3 global{options.nx, options.ny, options.nz};
        cfd::distributed::MpiCartesianRuntime runtime(global, 3, {true, true, true});

        int result = 0;
        if (options.backend == "cpu") {
            if (options.lattice == "d3q19") result = run_cpu<cfd::lbm::D3Q19Descriptor>(options, runtime, global, "d3q19");
            else result = run_cpu<cfd::lbm::D3Q27Descriptor>(options, runtime, global, "d3q27");
        } else {
#if defined(CFD_HAS_SYCL)
            if (options.lattice == "d3q19") result = run_sycl<cfd::lbm::D3Q19Descriptor>(options, runtime, global, "d3q19");
            else result = run_sycl<cfd::lbm::D3Q27Descriptor>(options, runtime, global, "d3q27");
#else
            throw std::runtime_error("SYCL backend requested but CFD_ENABLE_SYCL=OFF");
#endif
        }

        if (!options.csv) {
            std::cout << "rank=" << runtime.rank()
                      << " local_rank=" << runtime.local_rank() << '/' << runtime.local_size()
                      << " node=" << runtime.processor_name()
                      << " brick=" << runtime.brick().extent.x << 'x' << runtime.brick().extent.y << 'x' << runtime.brick().extent.z
                      << " begin=" << runtime.brick().begin.x << ',' << runtime.brick().begin.y << ',' << runtime.brick().begin.z
                      << '\n';
        }
        return result;
    } catch (const std::exception& e) {
        std::cerr << "cfd-distributed: " << e.what() << '\n';
        return 1;
    }
}

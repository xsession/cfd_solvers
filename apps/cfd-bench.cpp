#include "cfd/core/parallel.hpp"
#include "cfd/solvers/lbm/esoteric_pull.hpp"
#include "cfd/solvers/lbm/one_step_pull.hpp"
#if defined(CFD_HAS_SYCL)
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#endif

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

struct Options {
    std::string lattice{"d3q19"};
    std::string backend{"cpu"};
    std::string streaming{"in-place"};
    std::size_t nx{64};
    std::size_t ny{64};
    std::size_t nz{64};
    std::size_t warmup{10};
    std::size_t steps{100};
    unsigned threads{0};
    float tau{0.6F};
    bool csv{false};
};

std::size_t parse_size(const char* text, std::string_view option) {
    const auto value = std::strtoull(text, nullptr, 10);
    if (value == 0ULL) throw std::invalid_argument(std::string(option) + " must be > 0");
    return static_cast<std::size_t>(value);
}

float parse_float(const char* text, std::string_view option) {
    const float value = std::strtof(text, nullptr);
    if (!(value > 0.5F)) throw std::invalid_argument(std::string(option) + " must be > 0.5");
    return value;
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        const auto require_value = [&](std::string_view name) -> const char* {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(name) + " requires a value");
            return argv[++i];
        };
        if (arg == "--lattice") options.lattice = require_value(arg);
        else if (arg == "--backend") options.backend = require_value(arg);
        else if (arg == "--streaming") options.streaming = require_value(arg);
        else if (arg == "--nx") options.nx = parse_size(require_value(arg), arg);
        else if (arg == "--ny") options.ny = parse_size(require_value(arg), arg);
        else if (arg == "--nz") options.nz = parse_size(require_value(arg), arg);
        else if (arg == "--warmup") options.warmup = parse_size(require_value(arg), arg);
        else if (arg == "--steps") options.steps = parse_size(require_value(arg), arg);
        else if (arg == "--threads") options.threads = static_cast<unsigned>(parse_size(require_value(arg), arg));
        else if (arg == "--tau") options.tau = parse_float(require_value(arg), arg);
        else if (arg == "--csv") options.csv = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout
                << "Usage: cfd-bench [options]\n"
                << "  --lattice d2q9|d3q19|d3q27\n"
                << "  --backend cpu|sycl\n"
                << "  --streaming in-place|pull|both   (pull/both are CPU-only)\n"
                << "  --nx N --ny N --nz N\n"
                << "  --warmup N --steps N --tau T\n"
                << "  --threads N      OpenMP thread count (CPU)\n"
                << "  --csv            CSV output\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + std::string(arg));
        }
    }
    return options;
}

template<class Descriptor>
void print_result(const Options& options,
                  std::string_view backend,
                  std::string_view streaming,
                  std::string_view device,
                  double seconds,
                  std::size_t cells,
                  std::size_t allocated_population_bytes,
                  double mass_error_abs,
                  double mass_reference,
                  bool print_header = true) {
    const double updates = static_cast<double>(cells) * static_cast<double>(options.steps);
    const double mlups = updates / seconds / 1.0e6;

    // Each BGK step must at minimum read and write every population once. This is
    // a useful algorithmic traffic estimate; actual hardware traffic may be higher.
    const double bytes_per_update = 2.0 * static_cast<double>(Descriptor::q) * sizeof(float);
    const double estimated_gbps = mlups * bytes_per_update / 1000.0;
    const double allocated_mib = static_cast<double>(allocated_population_bytes) / (1024.0 * 1024.0);
    const double single_grid_mib = static_cast<double>(Descriptor::q) * static_cast<double>(cells) *
                                   sizeof(float) / (1024.0 * 1024.0);
    const double mass_error_rel = mass_reference != 0.0 ? mass_error_abs / std::abs(mass_reference) : 0.0;

    if (options.csv) {
        if (print_header) {
            std::cout << "lattice,backend,streaming,device,nx,ny,nz,q,steps,seconds,mlups,"
                         "estimated_ddf_gbps,allocated_ddf_mib,single_grid_ddf_mib,mass_error_abs,mass_error_rel\n";
        }
        std::cout << options.lattice << ',' << backend << ',' << streaming << ",\"" << device << "\"," << options.nx << ','
                  << options.ny << ',' << options.nz << ',' << Descriptor::q << ',' << options.steps << ','
                  << std::setprecision(9) << seconds << ',' << mlups << ',' << estimated_gbps << ','
                  << allocated_mib << ',' << single_grid_mib << ',' << mass_error_abs << ',' << mass_error_rel << '\n';
        return;
    }

    std::cout << "lattice=" << options.lattice
              << " backend=" << backend
              << " streaming=" << streaming
              << " device=\"" << device << "\"\n"
              << "grid=" << options.nx << 'x' << options.ny << 'x' << options.nz
              << " q=" << Descriptor::q << " steps=" << options.steps << '\n'
              << std::fixed << std::setprecision(3)
              << "time_s=" << seconds
              << " MLUPS=" << mlups
              << " estimated_DDF_GB/s=" << estimated_gbps << '\n'
              << "allocated_DDF_MiB=" << allocated_mib
              << " single_grid_DDF_MiB=" << single_grid_mib << '\n'
              << std::scientific << "mass_error_abs=" << mass_error_abs
              << " mass_error_rel=" << mass_error_rel << '\n';
}

template<class Descriptor>
int run_cpu_inplace(const Options& options, bool print_header = true) {
    using Solver = cfd::lbm::EsotericPullSolver<Descriptor>;
    if (options.threads > 0U) cfd::core::set_cpu_threads(options.threads);
    const cfd::lbm::InPlaceLbmConfig config{options.nx, options.ny, options.nz, options.tau};
    Solver solver(config);
    solver.initialize_taylor_green(0.02F);
    const double mass0 = solver.mass();
    solver.step(options.warmup);

    const auto begin = std::chrono::steady_clock::now();
    solver.step(options.steps);
    const auto end = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(end - begin).count();
    const double mass_error = std::abs(solver.mass() - mass0);
    print_result<Descriptor>(options, "cpu", "in-place", "OpenMP/host", seconds, solver.cells(),
                             solver.population_bytes(), mass_error, mass0, print_header);
    return 0;
}

template<class Descriptor>
int run_cpu_pull(const Options& options, bool print_header = true) {
    using Solver = cfd::lbm::OneStepPullSolver<Descriptor>;
    if (options.threads > 0U) cfd::core::set_cpu_threads(options.threads);
    const cfd::lbm::InPlaceLbmConfig config{options.nx, options.ny, options.nz, options.tau};
    Solver solver(config);
    solver.initialize_taylor_green(0.02F);
    const double mass0 = solver.mass();
    solver.step(options.warmup);

    const auto begin = std::chrono::steady_clock::now();
    solver.step(options.steps);
    const auto end = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(end - begin).count();
    const double mass_error = std::abs(solver.mass() - mass0);
    print_result<Descriptor>(options, "cpu", "two-grid-pull", "OpenMP/host", seconds, solver.cells(),
                             solver.population_bytes(), mass_error, mass0, print_header);
    return 0;
}

#if defined(CFD_HAS_SYCL)
template<class Descriptor>
int run_sycl(const Options& options) {
    using Solver = cfd::lbm::EsotericPullSyclSolver<Descriptor>;
    const cfd::lbm::InPlaceLbmConfig config{options.nx, options.ny, options.nz, options.tau};
    Solver solver(config);
    solver.initialize_taylor_green(0.02F);
    const auto initial = solver.download_macroscopic();
    double mass0 = 0.0;
    for (float rho : initial.rho) mass0 += rho;
    solver.step(options.warmup);
    solver.wait();

    const auto begin = std::chrono::steady_clock::now();
    solver.step(options.steps);
    solver.wait();
    const auto end = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(end - begin).count();
    const auto final = solver.download_macroscopic();
    double mass1 = 0.0;
    for (float rho : final.rho) mass1 += rho;
    print_result<Descriptor>(options, "sycl", "in-place", solver.device_name(), seconds, solver.cells(),
                             solver.population_bytes(), std::abs(mass1 - mass0), mass0);
    return 0;
}
#endif

template<class Descriptor>
int dispatch_backend(const Options& options) {
    if (options.backend == "cpu") {
        if (options.streaming == "in-place") return run_cpu_inplace<Descriptor>(options);
        if (options.streaming == "pull") return run_cpu_pull<Descriptor>(options);
        if (options.streaming == "both") {
            run_cpu_inplace<Descriptor>(options, true);
            return run_cpu_pull<Descriptor>(options, !options.csv);
        }
        throw std::invalid_argument("streaming must be in-place, pull, or both");
    }
#if defined(CFD_HAS_SYCL)
    if (options.backend == "sycl") {
        if (options.streaming != "in-place") {
            throw std::invalid_argument("SYCL benchmark currently supports only --streaming in-place");
        }
        return run_sycl<Descriptor>(options);
    }
#else
    if (options.backend == "sycl") {
        std::cerr << "SYCL backend was not compiled; configure with -DCFD_ENABLE_SYCL=ON\n";
        return 2;
    }
#endif
    throw std::invalid_argument("backend must be cpu or sycl");
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options = parse_options(argc, argv);
        if (options.lattice == "d2q9") {
            options.nz = 1;
            return dispatch_backend<cfd::lbm::D2Q9InPlaceDescriptor>(options);
        }
        if (options.lattice == "d3q19") return dispatch_backend<cfd::lbm::D3Q19Descriptor>(options);
        if (options.lattice == "d3q27") return dispatch_backend<cfd::lbm::D3Q27Descriptor>(options);
        throw std::invalid_argument("lattice must be d2q9, d3q19, or d3q27");
    } catch (const std::exception& error) {
        std::cerr << "cfd-bench: " << error.what() << '\n';
        return 1;
    }
}

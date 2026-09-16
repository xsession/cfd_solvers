#include "cfd/distributed/selective_halo.hpp"
#include "cfd/solvers/lbm/descriptors.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
struct Options {
    std::size_t nx{64};
    std::size_t ny{64};
    std::size_t nz{64};
    std::string lattice{"d3q19"};
};

std::size_t parse_size(std::string_view value) {
    const auto n = static_cast<std::size_t>(std::stoull(std::string(value)));
    if (n == 0) throw std::invalid_argument("extent must be positive");
    return n;
}

Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        auto value = [&]() -> std::string_view {
            if (i + 1 >= argc) throw std::invalid_argument(std::string(arg) + " requires a value");
            return argv[++i];
        };
        if (arg == "--nx") o.nx = parse_size(value());
        else if (arg == "--ny") o.ny = parse_size(value());
        else if (arg == "--nz") o.nz = parse_size(value());
        else if (arg == "--lattice") o.lattice = std::string(value());
        else if (arg == "--help") {
            std::cout << "Usage: cfd-halo-plan [--lattice d3q19|d3q27 --nx N --ny N --nz N]\n";
            std::exit(0);
        } else throw std::invalid_argument("unknown option: " + std::string(arg));
    }
    if (o.lattice != "d3q19" && o.lattice != "d3q27") throw std::invalid_argument("invalid lattice");
    return o;
}

template<class Descriptor>
void report(std::string_view name, cfd::distributed::Extent3 e) {
    const auto full = cfd::distributed::full_halo_traffic<Descriptor>(e);
    const auto compact = cfd::distributed::selective_halo_traffic<Descriptor>(e);
    const double ratio = full.bytes == 0 ? 0.0 : static_cast<double>(compact.bytes) / static_cast<double>(full.bytes);
    std::cout << "lattice=" << name
              << " brick=" << e.x << 'x' << e.y << 'x' << e.z
              << " full_messages=" << full.messages
              << " selective_messages=" << compact.messages
              << " full_bytes=" << full.bytes
              << " selective_bytes=" << compact.bytes
              << " traffic_ratio=" << std::fixed << std::setprecision(6) << ratio
              << " reduction_percent=" << std::setprecision(2) << (100.0 * (1.0 - ratio)) << '\n';
}
} // namespace

int main(int argc, char** argv) {
    try {
        const auto o = parse_options(argc, argv);
        const cfd::distributed::Extent3 e{o.nx, o.ny, o.nz};
        if (o.lattice == "d3q19") report<cfd::lbm::D3Q19Descriptor>(o.lattice, e);
        else report<cfd::lbm::D3Q27Descriptor>(o.lattice, e);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "cfd-halo-plan: " << e.what() << '\n';
        return 1;
    }
}

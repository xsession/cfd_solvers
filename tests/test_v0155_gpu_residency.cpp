#include "cfd/core/device_residency.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/sycl_sparse.hpp"
#include "cfd/solvers/lbm/esoteric_pull_sycl.hpp"
#endif

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition) throw std::runtime_error(std::string(message));
}

void test_transfer_ledger() {
    cfd::core::DeviceTransferStats stats;
    require(stats.host_transfer_bytes() == 0U, "new transfer ledger is empty");
    stats.record_host_to_device(128U);
    stats.record_device_to_host(64U);
    stats.record_device_to_device(32U);
    stats.record_synchronization(2U);
    require(stats.host_transfer_bytes() == 192U, "host traffic is counted separately");
    require(stats.explicit_transfer_bytes() == 224U, "all explicit traffic is counted");
    require(stats.synchronization_points == 2U, "synchronization points are counted");
    stats.reset();
    require(stats.explicit_transfer_bytes() == 0U && stats.synchronization_points == 0U,
            "transfer ledger reset");
}

#if defined(CFD_HAS_SYCL)

void test_sparse_direct_device_seam() {
    cfd::core::CsrBuilder builder(3U, 3U);
    builder.add(0U, 0U, 4.0); builder.add(0U, 1U, -1.0);
    builder.add(1U, 0U, -1.0); builder.add(1U, 1U, 4.0); builder.add(1U, 2U, -1.0);
    builder.add(2U, 1U, -1.0); builder.add(2U, 2U, 3.0);
    const auto matrix = builder.build();
    cfd::core::SyclCsrLinearAlgebra solver(matrix);

    const std::vector<double> exact{1.0, 2.0, 3.0};
    std::vector<double> rhs(3U, 0.0);
    matrix.multiply(exact, rhs);

    auto& queue = solver.queue();
    double* d_rhs = sycl::malloc_device<double>(3U, queue);
    double* d_x = sycl::malloc_device<double>(3U, queue);
    if (!d_rhs || !d_x) {
        if (d_rhs) sycl::free(d_rhs, queue);
        if (d_x) sycl::free(d_x, queue);
        throw std::bad_alloc{};
    }
    try {
        const std::vector<double> zero(3U, 0.0);
        queue.memcpy(d_rhs, rhs.data(), 3U * sizeof(double));
        queue.memcpy(d_x, zero.data(), 3U * sizeof(double)).wait_and_throw();
        solver.reset_transfer_stats();

        const auto result = solver.conjugate_gradient_device(d_rhs, d_x, 64U, 1.0e-12);
        require(result.converged, "direct-device SYCL CG converges");
        require(solver.transfer_stats().host_transfer_bytes() == 0U,
                "direct-device SYCL CG does not stage volume vectors through host memory");

        std::vector<double> x(3U, 0.0);
        queue.memcpy(x.data(), d_x, 3U * sizeof(double)).wait_and_throw();
        for (std::size_t i = 0U; i < x.size(); ++i) {
            require(std::abs(x[i] - exact[i]) < 1.0e-9,
                    "direct-device SYCL CG solution parity");
        }
    } catch (...) {
        sycl::free(d_rhs, queue);
        sycl::free(d_x, queue);
        throw;
    }
    sycl::free(d_rhs, queue);
    sycl::free(d_x, queue);
}

void test_lbm_device_residency() {
    constexpr std::size_t nx = 12U;
    constexpr std::size_t ny = 10U;
    constexpr std::size_t nz = 8U;
    cfd::lbm::D3Q19SyclSolver solver({nx, ny, nz, 0.72F});
    solver.initialize_taylor_green(0.01F);
    solver.wait();
    solver.reset_transfer_stats();

    const double mass0 = solver.total_mass();
    solver.step(5U);
    solver.wait();
    const double mass1 = solver.total_mass();
    require(solver.transfer_stats().host_transfer_bytes() == 0U,
            "device-resident LBM stepping and scalar reduction do not transfer the lattice");
    require(mass0 > 0.0 && mass1 > 0.0, "device mass reduction returns finite positive mass");

    const auto fields = solver.download_macroscopic();
    require(fields.rho.size() == nx * ny * nz, "macroscopic extraction size");
    require(solver.transfer_stats().device_to_host_bytes ==
                4U * nx * ny * nz * sizeof(float),
            "macroscopic extraction transfers only rho/u, not all populations");
}
#endif

} // namespace

int main() {
    try {
        test_transfer_ledger();
#if defined(CFD_HAS_SYCL)
        test_sparse_direct_device_seam();
        test_lbm_device_residency();
#endif
        std::cout << "v0.15.5 GPU residency tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.5 GPU residency test failure: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/numa.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/core/sycl_sparse.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

#if defined(CFD_HAS_SYCL)
cfd::core::CsrMatrix make_spd_matrix() {
    cfd::core::CsrBuilder builder(3U, 3U);
    builder.add(0U, 0U, 4.0); builder.add(0U, 1U, -1.0);
    builder.add(1U, 0U, -1.0); builder.add(1U, 1U, 4.0); builder.add(1U, 2U, -1.0);
    builder.add(2U, 1U, -1.0); builder.add(2U, 2U, 3.0);
    return builder.build();
}
#endif

void test_numa_planning() {
    const auto parsed = cfd::core::parse_cpu_list("0-3,8,10-11\n");
    require(parsed == std::vector<unsigned>({0U, 1U, 2U, 3U, 8U, 10U, 11U}), "CPU-list parser regression");

    const auto nodes = cfd::core::discover_numa_nodes();
    require(!nodes.empty(), "NUMA topology should expose at least a fallback node");
    std::size_t available = 0U;
    for (const auto& node : nodes) {
        require(!node.cpus.empty(), "NUMA node must contain an allowed CPU");
        available += node.cpus.size();
    }
    require(available > 0U, "NUMA topology must expose CPUs");

    const std::size_t threads = std::min<std::size_t>(available, 4U);
    const auto compact = cfd::core::plan_numa_thread_placement(threads, cfd::core::NumaPlacementPolicy::compact);
    const auto spread = cfd::core::plan_numa_thread_placement(threads, cfd::core::NumaPlacementPolicy::spread);
    require(compact.thread_count() == threads && compact.thread_nodes.size() == threads, "compact placement size mismatch");
    require(spread.thread_count() == threads && spread.thread_nodes.size() == threads, "spread placement size mismatch");

    const auto explicit_node = cfd::core::plan_numa_thread_placement(
        threads, cfd::core::NumaPlacementPolicy::explicit_node, nodes.front().id);
    require(std::all_of(explicit_node.thread_nodes.begin(), explicit_node.thread_nodes.end(),
                        [&](int node) { return node == nodes.front().id; }),
            "explicit NUMA node plan escaped requested node");

    std::vector<double> touched(4096U, -1.0);
    cfd::core::numa_first_touch<double>(touched, 2.5);
    require(touched.front() == 2.5 && touched.back() == 2.5, "portable NUMA first-touch regression");
}

#if defined(CFD_HAS_SYCL)
void test_sycl_sparse_krylov() {
    const auto matrix = make_spd_matrix();
    cfd::core::SyclCsrLinearAlgebra sycl_matrix(matrix);
    const std::vector<double> exact{1.0, 2.0, 3.0};
    std::vector<double> rhs(3U);
    matrix.multiply(exact, rhs);

    std::vector<double> sycl_rhs(3U);
    sycl_matrix.multiply(exact, sycl_rhs);
    for (std::size_t i = 0U; i < rhs.size(); ++i) {
        require(std::abs(rhs[i] - sycl_rhs[i]) < 1.0e-12, "SYCL SpMV parity regression");
    }

    std::vector<double> x(3U, 0.0);
    const auto result = sycl_matrix.conjugate_gradient(rhs, x, 64U, 1.0e-12);
    require(result.converged, "SYCL CG failed to converge");
    for (std::size_t i = 0U; i < x.size(); ++i) {
        require(std::abs(x[i] - exact[i]) < 1.0e-9, "SYCL CG solution regression");
    }
}
#endif

} // namespace

int main() {
    try {
        test_numa_planning();
#if defined(CFD_HAS_SYCL)
        test_sycl_sparse_krylov();
#endif
        std::cout << "v0.14.0 common HPC runtime tests passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "v0.14.0 common HPC runtime test failure: " << e.what() << '\n';
        return 1;
    }
}

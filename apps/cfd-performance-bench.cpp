#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/distributed_csr.hpp"
#include "cfd/core/parallel.hpp"
#include "cfd/core/structured_stencil.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    std::size_t nx{32U};
    std::size_t ny{32U};
    std::size_t nz{32U};
    std::size_t warmup{5U};
    std::size_t steps{30U};
    unsigned threads{};
    bool quick{};
};

std::size_t parse_size(const char* text, std::string_view option) {
    const auto value = std::strtoull(text, nullptr, 10);
    if (value == 0ULL)
        throw std::invalid_argument(std::string(option) + " must be > 0");
    return static_cast<std::size_t>(value);
}

Options parse_options(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        const auto require_value = [&](std::string_view name) -> const char* {
            if (i + 1 >= argc)
                throw std::invalid_argument(std::string(name) + " requires a value");
            return argv[++i];
        };
        if (arg == "--nx")
            options.nx = parse_size(require_value(arg), arg);
        else if (arg == "--ny")
            options.ny = parse_size(require_value(arg), arg);
        else if (arg == "--nz")
            options.nz = parse_size(require_value(arg), arg);
        else if (arg == "--warmup")
            options.warmup = parse_size(require_value(arg), arg);
        else if (arg == "--steps")
            options.steps = parse_size(require_value(arg), arg);
        else if (arg == "--threads")
            options.threads = static_cast<unsigned>(parse_size(require_value(arg), arg));
        else if (arg == "--quick")
            options.quick = true;
        else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: cfd-performance-bench [options]\n"
                      << "  --nx N --ny N --nz N\n"
                      << "  --warmup N --steps N --threads N\n"
                      << "  --quick\n";
            std::exit(0);
        } else {
            throw std::invalid_argument("unknown option: " + std::string(arg));
        }
    }
    if (options.quick) {
        options.warmup = std::min(options.warmup, std::size_t{2U});
        options.steps = std::min(options.steps, std::size_t{8U});
    }
    return options;
}

std::size_t cell_index(std::size_t nx, std::size_t ny, std::size_t i, std::size_t j, std::size_t k) {
    return (k * ny + j) * nx + i;
}

std::size_t structured_nonzeros(std::size_t nx, std::size_t ny, std::size_t nz) {
    std::size_t result = 0U;
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                result += 1U + static_cast<std::size_t>(i > 0U) + static_cast<std::size_t>(i + 1U < nx) +
                          static_cast<std::size_t>(j > 0U) + static_cast<std::size_t>(j + 1U < ny) +
                          static_cast<std::size_t>(k > 0U) + static_cast<std::size_t>(k + 1U < nz);
            }
        }
    }
    return result;
}

cfd::core::CsrMatrix make_seven_point_matrix(std::size_t nx, std::size_t ny, std::size_t nz) {
    const std::size_t cells = nx * ny * nz;
    cfd::core::CsrBuilder builder(cells, cells);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto row = cell_index(nx, ny, i, j, k);
                builder.add(row, row, 6.0);
                if (i > 0U)
                    builder.add(row, cell_index(nx, ny, i - 1U, j, k), -1.0);
                if (i + 1U < nx)
                    builder.add(row, cell_index(nx, ny, i + 1U, j, k), -1.0);
                if (j > 0U)
                    builder.add(row, cell_index(nx, ny, i, j - 1U, k), -1.0);
                if (j + 1U < ny)
                    builder.add(row, cell_index(nx, ny, i, j + 1U, k), -1.0);
                if (k > 0U)
                    builder.add(row, cell_index(nx, ny, i, j, k - 1U), -1.0);
                if (k + 1U < nz)
                    builder.add(row, cell_index(nx, ny, i, j, k + 1U), -1.0);
            }
        }
    }
    return builder.build();
}

void fill_vector(std::vector<double>& values, double phase) {
    for (std::size_t i = 0; i < values.size(); ++i) {
        values[i] = 0.25 + 0.001 * static_cast<double>((i * 17U + static_cast<std::size_t>(phase * 100.0)) % 997U);
    }
}

double checksum(std::span<const double> values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

double max_difference(std::span<const double> a, std::span<const double> b) {
    if (a.size() != b.size())
        throw std::invalid_argument("benchmark comparison size mismatch");
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
        result = std::max(result, std::abs(a[i] - b[i]));
    return result;
}

double residual_l2(std::span<const double> rhs, std::span<const double> value) {
    if (rhs.size() != value.size())
        throw std::invalid_argument("benchmark residual size mismatch");
    double sum = 0.0;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const double difference = rhs[i] - value[i];
        sum += difference * difference;
    }
    return std::sqrt(sum);
}

std::size_t resident_bytes() {
#if defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.rfind("VmRSS:", 0U) == 0U) {
            std::istringstream fields(line.substr(6U));
            std::size_t value{};
            fields >> value;
            return value * 1024U;
        }
    }
#endif
    return 0U;
}

void legacy_distributed_multiply(const cfd::core::CsrMatrix& matrix, std::size_t row_begin, std::size_t row_end,
                                 std::span<const std::size_t> halo_columns, std::span<const double> local,
                                 std::span<const double> halo, std::span<double> y) {
    const auto& offsets = matrix.row_offsets();
    const auto& columns = matrix.column_indices();
    const auto& values = matrix.values();
    for (std::size_t local_row = 0; local_row < row_end - row_begin; ++local_row) {
        const std::size_t row = row_begin + local_row;
        double sum = 0.0;
        for (std::size_t k = offsets[row]; k < offsets[row + 1U]; ++k) {
            const auto column = columns[k];
            double x = 0.0;
            if (column >= row_begin && column < row_end) {
                x = local[column - row_begin];
            } else {
                const auto it = std::lower_bound(halo_columns.begin(), halo_columns.end(), column);
                if (it == halo_columns.end() || *it != column)
                    throw std::runtime_error("legacy halo lookup failed");
                x = halo[static_cast<std::size_t>(it - halo_columns.begin())];
            }
            sum += values[k] * x;
        }
        y[local_row] = sum;
    }
}

struct Timing {
    double elapsed_ms{};
    double checksum{};
};

template <class Function> Timing time_kernel(const Options& options, std::span<double> output, Function&& function) {
    for (std::size_t i = 0; i < options.warmup; ++i)
        function();
    const auto start = Clock::now();
    for (std::size_t i = 0; i < options.steps; ++i)
        function();
    const auto end = Clock::now();
    return {std::chrono::duration<double, std::milli>(end - start).count(), checksum(output)};
}

void emit_record(std::string_view benchmark, std::string_view implementation, const Options& options, std::size_t cells,
                 std::size_t nonzeros, std::size_t bytes_per_step, std::size_t working_set, const Timing& timing,
                 double oracle_error, double residual) {
    const double seconds = timing.elapsed_ms * 1.0e-3;
    const double work = static_cast<double>(cells) * static_cast<double>(options.steps);
    const double throughput = seconds > 0.0 ? work / seconds : 0.0;
    const double gbps = seconds > 0.0
                            ? static_cast<double>(bytes_per_step) * static_cast<double>(options.steps) / seconds / 1.0e9
                            : 0.0;
    const double flops =
        seconds > 0.0 ? 2.0 * static_cast<double>(nonzeros) * static_cast<double>(options.steps) / seconds : 0.0;
    const unsigned effective_threads = options.threads > 0U ? options.threads : cfd::core::cpu_thread_capacity();
    std::cout << std::setprecision(17) << "CFD_PERF {\"schema\":1" << ",\"benchmark\":\"" << benchmark << "\""
              << ",\"implementation\":\"" << implementation << "\"" << ",\"backend\":\"cpu\""
              << ",\"nx\":" << options.nx << ",\"ny\":" << options.ny << ",\"nz\":" << options.nz
              << ",\"cells\":" << cells << ",\"nonzeros\":" << nonzeros << ",\"threads\":" << effective_threads
              << ",\"warmup\":" << options.warmup << ",\"steps\":" << options.steps
              << ",\"elapsed_ms\":" << timing.elapsed_ms
              << ",\"ms_per_step\":" << timing.elapsed_ms / static_cast<double>(options.steps)
              << ",\"cells_per_second\":" << throughput << ",\"estimated_gbps\":" << gbps
              << ",\"estimated_flops_per_second\":" << flops << ",\"bytes_per_step\":" << bytes_per_step
              << ",\"working_set_bytes\":" << working_set << ",\"checksum\":" << timing.checksum
              << ",\"oracle_max_abs_error\":" << oracle_error << ",\"residual_l2\":" << residual << "}\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parse_options(argc, argv);
        if (options.threads > 0U)
            cfd::core::set_cpu_threads(options.threads);

        const std::size_t cells = options.nx * options.ny * options.nz;
        const std::size_t stencil_nonzeros = structured_nonzeros(options.nx, options.ny, options.nz);
        std::vector<double> x(cells), rhs(cells), output(cells), residual(cells), reference(cells);
        fill_vector(x, 0.0);
        fill_vector(rhs, 1.0);

        const auto structured_setup_start = Clock::now();
        const cfd::core::StructuredSevenPointOperator stencil(options.nx, options.ny, options.nz);
        const auto structured_setup_end = Clock::now();
        (void)structured_setup_start;
        (void)structured_setup_end;

        stencil.apply(x, reference);
        const Timing structured_apply = time_kernel(options, output, [&] { stencil.apply(x, output); });
        emit_record("structured_spmv", "structured_seven_point", options, cells, stencil_nonzeros,
                    stencil.estimated_apply_bytes(), resident_bytes(), structured_apply,
                    max_difference(reference, output), residual_l2(rhs, output));

        const Timing structured_explicit_residual = time_kernel(options, residual, [&] {
            stencil.apply(x, output);
            for (std::size_t i = 0; i < cells; ++i)
                residual[i] = rhs[i] - output[i];
        });
        emit_record("structured_residual", "structured_apply_plus_residual", options, cells, stencil_nonzeros,
                    stencil.estimated_apply_bytes() + cells * sizeof(double) * 2U, resident_bytes(),
                    structured_explicit_residual, 0.0, residual_l2(rhs, output));

        const Timing structured_fused_residual =
            time_kernel(options, residual, [&] { (void)stencil.apply_residual_l2(rhs, x, residual); });
        emit_record("structured_residual", "structured_fused_apply_residual", options, cells, stencil_nonzeros,
                    stencil.estimated_fused_residual_bytes(), resident_bytes(), structured_fused_residual, 0.0,
                    stencil.apply_residual_l2(rhs, x, residual));

        const auto csr_setup_start = Clock::now();
        const auto matrix = make_seven_point_matrix(options.nx, options.ny, options.nz);
        const auto csr_setup_end = Clock::now();
        matrix.multiply(x, reference);
        const Timing csr_apply = time_kernel(options, output, [&] { matrix.multiply(x, output); });
        emit_record("csr_spmv", "csr", options, cells, matrix.nonzeros(),
                    matrix.nonzeros() * (sizeof(double) + sizeof(std::size_t)) + (cells + 1U) * sizeof(std::size_t),
                    resident_bytes(), csr_apply, max_difference(reference, output), residual_l2(rhs, output));

        const std::size_t row_begin = 0U;
        const std::size_t row_end = cells / 2U;
        cfd::core::DistributedCsrPartition partition(matrix, row_begin, row_end);
        std::vector<double> local_x(row_end - row_begin), halo_x(partition.halo_columns().size()),
            local_output(row_end - row_begin), legacy_output(row_end - row_begin);
        std::copy(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(row_end), local_x.begin());
        for (std::size_t i = 0; i < partition.halo_columns().size(); ++i)
            halo_x[i] = x[partition.halo_columns()[i]];
        legacy_distributed_multiply(matrix, row_begin, row_end, partition.halo_columns(), local_x, halo_x,
                                    legacy_output);
        partition.multiply(local_x, halo_x, local_output);
        const double distributed_error = max_difference(legacy_output, local_output);
        const Timing legacy = time_kernel(options, legacy_output, [&] {
            legacy_distributed_multiply(matrix, row_begin, row_end, partition.halo_columns(), local_x, halo_x,
                                        legacy_output);
        });
        emit_record("distributed_csr_spmv", "legacy_lower_bound", options, row_end - row_begin, matrix.nonzeros(),
                    matrix.nonzeros() * (sizeof(double) + sizeof(std::size_t)), resident_bytes(), legacy,
                    distributed_error, 0.0);

        const Timing remapped =
            time_kernel(options, local_output, [&] { partition.multiply(local_x, halo_x, local_output); });
        emit_record("distributed_csr_spmv", "remapped_parallel", options, row_end - row_begin, matrix.nonzeros(),
                    matrix.nonzeros() * (sizeof(double) + sizeof(std::size_t)), resident_bytes(), remapped,
                    distributed_error, 0.0);

        const double setup_structured_ms =
            std::chrono::duration<double, std::milli>(structured_setup_end - structured_setup_start).count();
        const double setup_csr_ms = std::chrono::duration<double, std::milli>(csr_setup_end - csr_setup_start).count();
        std::cout << std::setprecision(17) << "CFD_PERF_SETUP {\"schema\":1,\"nx\":" << options.nx
                  << ",\"ny\":" << options.ny << ",\"nz\":" << options.nz
                  << ",\"structured_setup_ms\":" << setup_structured_ms << ",\"csr_setup_ms\":" << setup_csr_ms
                  << ",\"csr_nonzeros\":" << matrix.nonzeros()
                  << ",\"distributed_halo_values\":" << partition.halo_columns().size() << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "cfd-performance-bench: " << error.what() << '\n';
        return 1;
    }
}

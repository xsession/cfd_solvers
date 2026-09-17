#include "cfd/fvm/resident_sycl.hpp"
#include "cfd/fvm/operators.hpp"

#if defined(CFD_HAS_SYCL)
#include "cfd/core/sycl_sparse.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

void require(bool condition, std::string_view message) {
    if (!condition) throw std::runtime_error(std::string(message));
}

double csr_value(const cfd::core::CsrMatrix& a, std::size_t row, std::size_t col) {
    for (std::size_t k = a.row_offsets()[row]; k < a.row_offsets()[row + 1U]; ++k) {
        if (a.column_indices()[k] == col) return a.values()[k];
    }
    return 0.0;
}

void test_pressure_matrix_cpu_baseline() {
    const auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(4U, 3U, 2U, 1.0, 0.8, 0.6, 0.17);
    const auto matrix = cfd::fvm::build_orthogonal_pressure_matrix(mesh, 0.75, true);
    require(matrix.rows() == mesh.cell_count() && matrix.cols() == mesh.cell_count(),
            "resident FVM pressure matrix dimensions");
    for (std::size_t r = 0U; r < matrix.rows(); ++r) {
        require(csr_value(matrix, r, r) > 0.0, "pressure matrix has positive diagonal");
        for (std::size_t k = matrix.row_offsets()[r]; k < matrix.row_offsets()[r + 1U]; ++k) {
            const std::size_t c = matrix.column_indices()[k];
            require(std::abs(matrix.values()[k] - csr_value(matrix, c, r)) < 1.0e-12,
                    "orthogonal pressure matrix is symmetric");
        }
    }
    std::vector<double> x(matrix.rows());
    for (std::size_t i = 0U; i < x.size(); ++i) x[i] = 0.1 + std::sin(0.31 * static_cast<double>(i));
    std::vector<double> y(x.size(), 0.0);
    matrix.multiply(x, y);
    double energy = 0.0;
    for (std::size_t i = 0U; i < x.size(); ++i) energy += x[i] * y[i];
    require(energy > 0.0, "Dirichlet orthogonal pressure matrix is positive definite");

    const auto pinned = cfd::fvm::build_orthogonal_pressure_matrix(mesh, 0.75, false, 0U);
    require(std::abs(csr_value(pinned, 0U, 0U) - 1.0) < 1.0e-14,
            "pinned Neumann pressure reference has unit diagonal");
    for (std::size_t c = 1U; c < pinned.cols(); ++c) {
        require(std::abs(csr_value(pinned, 0U, c)) < 1.0e-14 &&
                std::abs(csr_value(pinned, c, 0U)) < 1.0e-14,
                "pinned pressure reference removes row and column couplings symmetrically");
    }
    std::vector<double> xp(pinned.rows(), 0.0), yp(pinned.rows(), 0.0);
    for (std::size_t i = 1U; i < xp.size(); ++i) xp[i] = 0.2 + std::cos(0.27 * static_cast<double>(i));
    pinned.multiply(xp, yp);
    double pinned_energy = 0.0;
    for (std::size_t i = 0U; i < xp.size(); ++i) pinned_energy += xp[i] * yp[i];
    require(pinned_energy > 0.0, "pinned Neumann pressure matrix is positive definite");
}

#if defined(CFD_HAS_SYCL)

void require_close(double a, double b, double tolerance, std::string_view what) {
    if (std::abs(a - b) > tolerance * std::max({1.0, std::abs(a), std::abs(b)})) {
        throw std::runtime_error(std::string(what));
    }
}

void test_resident_fvm_operators_and_sparse_seam() {
    const auto mesh = cfd::fvm::make_sheared_cartesian_hexa_mesh(5U, 4U, 3U, 1.2, 0.9, 0.7, 0.13);
    cfd::fvm::ResidentPolyMeshSycl resident(mesh);
    auto& q = resident.queue();
    const std::size_t nc = mesh.cell_count(), nf = mesh.face_count();

    double* d_phi = sycl::malloc_device<double>(nc, q);
    double* d_grad = sycl::malloc_device<double>(3U * nc, q);
    double* d_lap = sycl::malloc_device<double>(nc, q);
    double* d_vec = sycl::malloc_device<double>(3U * nc, q);
    double* d_div = sycl::malloc_device<double>(nc, q);
    double* d_flux = sycl::malloc_device<double>(nf, q);
    double* d_flux_div = sycl::malloc_device<double>(nc, q);
    double* d_mobility = sycl::malloc_device<double>(nc, q);
    double* d_velocity = sycl::malloc_device<double>(3U * nc, q);
    require(d_phi && d_grad && d_lap && d_vec && d_div && d_flux && d_flux_div && d_mobility && d_velocity,
            "resident FVM device allocations succeed");

    auto cleanup = [&]() {
        if (d_velocity) sycl::free(d_velocity, q);
        if (d_mobility) sycl::free(d_mobility, q);
        if (d_flux_div) sycl::free(d_flux_div, q);
        if (d_flux) sycl::free(d_flux, q);
        if (d_div) sycl::free(d_div, q);
        if (d_vec) sycl::free(d_vec, q);
        if (d_lap) sycl::free(d_lap, q);
        if (d_grad) sycl::free(d_grad, q);
        if (d_phi) sycl::free(d_phi, q);
    };

    try {
        std::vector<double> phi(nc), mobility(nc, 0.7);
        std::vector<cfd::fvm::Vec3> vec(nc);
        for (std::size_t c = 0U; c < nc; ++c) {
            const auto p = mesh.cells()[c].center;
            phi[c] = 0.3 + 1.2 * p.x - 0.4 * p.y + 0.2 * p.z + 0.1 * p.x * p.x;
            vec[c] = {0.5 + 0.2 * p.x, -0.1 + 0.3 * p.y, 0.15 * p.z};
        }
        resident.upload_cell_scalar(phi, d_phi);
        resident.upload_cell_vector(vec, d_vec);
        resident.upload_cell_scalar(mobility, d_mobility);
        resident.reset_transfer_stats();

        resident.gauss_gradient_scalar(d_phi, d_grad);
        resident.orthogonal_laplacian_scalar(d_phi, 0.42, d_lap);
        resident.gauss_divergence_vector(d_vec, d_div);
        resident.predictor_face_flux(d_vec, d_flux);
        resident.divergence_face_flux(d_flux, d_flux_div);
        resident.correct_velocity_from_pressure_gradient(d_vec, d_mobility, d_grad, d_velocity);
        resident.wait();
        require(resident.transfer_stats().host_transfer_bytes() == 0U,
                "resident FVM operator chain performs no host field transfer");
        require(resident.resident_bytes() > (4U * nc + 6U * nf) * sizeof(double),
                "resident FVM accounts for connectivity and operator metrics");

        std::vector<cfd::fvm::Vec3> grad(nc), corrected_velocity(nc);
        std::vector<double> lap(nc), div(nc), flux_div(nc);
        resident.download_cell_vector(d_grad, grad);
        resident.download_cell_scalar(d_lap, lap);
        resident.download_cell_scalar(d_div, div);
        resident.download_cell_scalar(d_flux_div, flux_div);
        resident.download_cell_vector(d_velocity, corrected_velocity);

        const auto grad_ref = cfd::fvm::gauss_gradient_scalar(mesh, phi);
        const auto lap_ref = cfd::fvm::orthogonal_laplacian_scalar(mesh, phi, 0.42);
        const auto div_ref = cfd::fvm::gauss_divergence_vector(mesh, vec);
        for (std::size_t c = 0U; c < nc; ++c) {
            require_close(grad[c].x, grad_ref[c].x, 1.0e-11, "resident FVM gradient x parity");
            require_close(grad[c].y, grad_ref[c].y, 1.0e-11, "resident FVM gradient y parity");
            require_close(grad[c].z, grad_ref[c].z, 1.0e-11, "resident FVM gradient z parity");
            require_close(lap[c], lap_ref[c], 1.0e-11, "resident FVM Laplacian parity");
            require_close(div[c], div_ref[c], 1.0e-11, "resident FVM divergence parity");
            require_close(flux_div[c], div_ref[c], 1.0e-11, "resident face-flux divergence parity");
            require_close(corrected_velocity[c].x, vec[c].x - mobility[c] * grad_ref[c].x, 1.0e-11,
                          "resident pressure velocity correction x");
        }

        // Reuse exactly the same queue/context in the persistent sparse solver.
        // This is the critical seam that lets pressure vectors remain resident.
        const auto pressure_matrix = cfd::fvm::build_orthogonal_pressure_matrix(mesh, mobility, true);
        cfd::core::SyclCsrLinearAlgebra pressure_solver(pressure_matrix, q);
        double* d_rhs = sycl::malloc_device<double>(nc, q);
        double* d_p = sycl::malloc_device<double>(nc, q);
        require(d_rhs && d_p, "resident pressure vectors allocate");
        try {
            std::vector<double> exact(nc), rhs(nc), zero(nc, 0.0), solved(nc);
            for (std::size_t c = 0U; c < nc; ++c) exact[c] = 0.05 + std::cos(0.19 * static_cast<double>(c));
            pressure_matrix.multiply(exact, rhs);
            resident.upload_cell_scalar(rhs, d_rhs); resident.upload_cell_scalar(zero, d_p);
            pressure_solver.reset_transfer_stats();
            const auto result = pressure_solver.conjugate_gradient_device(d_rhs, d_p, 600U, 1.0e-12);
            require(result.converged, "resident FVM direct-device pressure CG converges");
            require(pressure_solver.transfer_stats().host_transfer_bytes() == 0U,
                    "resident FVM pressure solve stages no volume vector through host");
            resident.download_cell_scalar(d_p, solved);
            for (std::size_t c = 0U; c < nc; ++c) {
                require_close(solved[c], exact[c], 2.0e-9, "resident FVM pressure CG solution parity");
            }
        } catch (...) {
            if (d_rhs) sycl::free(d_rhs, q);
            if (d_p) sycl::free(d_p, q);
            throw;
        }
        sycl::free(d_rhs, q); sycl::free(d_p, q);

        cfd::fvm::ResidentPressureProjectionSycl projection(mesh, 0.7, 0U, q);
        projection.set_predictor(vec);
        projection.reset_transfer_stats();
        const auto projected = projection.project(600U, 1.0e-11);
        require(projected.converged, "resident pressure projection CG converges");
        require(projection.hot_loop_host_transfer_bytes() == 0U,
                "resident pressure projection performs no bulk host transfer");
        const double continuity = projection.continuity_l2();
        require(std::isfinite(continuity) && continuity < 1.0e-8,
                "resident pressure projection reduces continuity residual");
    } catch (...) {
        cleanup(); throw;
    }
    cleanup();
}
#endif

} // namespace

int main() {
    try {
        test_pressure_matrix_cpu_baseline();
#if defined(CFD_HAS_SYCL)
        test_resident_fvm_operators_and_sparse_seam();
#endif
        std::cout << "v0.15.8 GPU FVM residency tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "v0.15.8 GPU FVM residency test failure: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

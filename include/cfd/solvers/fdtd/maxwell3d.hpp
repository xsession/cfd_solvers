#pragma once

#include "cfd/core/aligned_allocator.hpp"

#include <cstddef>

namespace cfd::fdtd {

enum class Boundary3D : unsigned char { pec, pmc };

struct Maxwell3DConfig {
    std::size_t nx{32};
    std::size_t ny{32};
    std::size_t nz{32};
    double dx{1.0e-3};
    double dy{1.0e-3};
    double dz{1.0e-3};
    double courant{0.95};
    double epsilon_r{1.0};
    double mu_r{1.0};
    Boundary3D boundary{Boundary3D::pec};
};

// 3-D Cartesian Yee scheme. The component staggering is implicit in the
// forward-H/backward-E difference pairing; arrays use a common padded logical
// extent to keep indexing and future accelerator kernels compact.
class Maxwell3D {
public:
    explicit Maxwell3D(Maxwell3DConfig config = {});

    void initialize_gaussian_ez(double amplitude = 1.0, double width_fraction = 0.12);
    void add_soft_ez_source(std::size_t i, std::size_t j, std::size_t k, double value);
    void step(std::size_t count = 1U);

    [[nodiscard]] double energy() const;
    [[nodiscard]] double max_field() const;
    [[nodiscard]] double dt() const noexcept { return dt_; }
    [[nodiscard]] std::size_t cell_count() const noexcept { return config_.nx*config_.ny*config_.nz; }
    [[nodiscard]] std::size_t nx() const noexcept { return config_.nx; }
    [[nodiscard]] std::size_t ny() const noexcept { return config_.ny; }
    [[nodiscard]] std::size_t nz() const noexcept { return config_.nz; }
    [[nodiscard]] double dx() const noexcept { return config_.dx; }
    [[nodiscard]] double dy() const noexcept { return config_.dy; }
    [[nodiscard]] double dz() const noexcept { return config_.dz; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& ex() const noexcept { return ex_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& ey() const noexcept { return ey_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& ez() const noexcept { return ez_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& hx() const noexcept { return hx_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& hy() const noexcept { return hy_; }
    [[nodiscard]] const cfd::core::AlignedVector<double>& hz() const noexcept { return hz_; }

private:
    Maxwell3DConfig config_;
    double dt_{};
    double epsilon_{};
    double mu_{};
    cfd::core::AlignedVector<double> ex_,ey_,ez_,hx_,hy_,hz_;

    [[nodiscard]] std::size_t idx(std::size_t i,std::size_t j,std::size_t k) const noexcept;
    void enforce_boundary();
};

} // namespace cfd::fdtd

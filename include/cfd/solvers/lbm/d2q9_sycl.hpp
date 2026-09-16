#pragma once

#include "cfd/solvers/lbm/d2q9.hpp"

#if defined(CFD_HAS_SYCL)
#include <sycl/sycl.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cfd::lbm {

class D2Q9SyclSolver {
public:
    explicit D2Q9SyclSolver(D2Q9Config config);
    ~D2Q9SyclSolver();

    D2Q9SyclSolver(const D2Q9SyclSolver&) = delete;
    D2Q9SyclSolver& operator=(const D2Q9SyclSolver&) = delete;

    void initialize_uniform(float rho = 1.0F, float ux = 0.0F, float uy = 0.0F);
    void initialize_taylor_green(float amplitude = 0.03F);
    void step(std::size_t count = 1);
    void wait();

    [[nodiscard]] std::vector<float> download_density() const;
    [[nodiscard]] std::string device_name() const;

private:
    D2Q9Config config_;
    std::size_t cells_{};
    sycl::queue queue_;
    float* f_{nullptr};
    float* next_{nullptr};
    std::uint8_t* solid_{nullptr};

    void upload_equilibrium(const std::vector<float>& rho,
                            const std::vector<float>& ux,
                            const std::vector<float>& uy);
};

} // namespace cfd::lbm
#endif

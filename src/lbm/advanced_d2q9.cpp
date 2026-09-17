#include "cfd/solvers/lbm/advanced_d2q9.hpp"
#include "cfd/core/parallel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::lbm {
namespace {

template<class Initializer>
void initialize_populations(const cfd::core::Grid2D& grid,
                            cfd::core::StaticSoA<float, 9>& f,
                            Initializer&& initializer) {
    cfd::core::parallel_for(grid.cells(), [&](std::size_t n) {
        const std::size_t x = n % grid.nx;
        const std::size_t y = n / grid.nx;
        const auto state = initializer(x, y);
        for (int d = 0; d < 9; ++d) {
            f(static_cast<std::size_t>(d), n) =
                D2Q9Descriptor::equilibrium(d, state[0], state[1], state[2]);
        }
    });
}

std::array<float, 3> macroscopic(const std::array<float, 9>& f) {
    float rho = 0.0F;
    float ux = 0.0F;
    float uy = 0.0F;
    for (int d = 0; d < 9; ++d) {
        const float value = f[static_cast<std::size_t>(d)];
        rho += value;
        ux += value * static_cast<float>(D2Q9Descriptor::cx(d));
        uy += value * static_cast<float>(D2Q9Descriptor::cy(d));
    }
    if (rho > 0.0F) {
        ux /= rho;
        uy /= rho;
    }
    return {rho, ux, uy};
}

std::size_t periodic_source(const cfd::core::Grid2D& grid,
                            std::size_t x,
                            std::size_t y,
                            int direction) noexcept {
    const auto nx = static_cast<long long>(grid.nx);
    const auto ny = static_cast<long long>(grid.ny);
    const auto sx = (static_cast<long long>(x) - D2Q9Descriptor::cx(direction) + nx) % nx;
    const auto sy = (static_cast<long long>(y) - D2Q9Descriptor::cy(direction) + ny) % ny;
    return static_cast<std::size_t>(sy) * grid.nx + static_cast<std::size_t>(sx);
}

template<class Field>
double mass_impl(const Field& f, std::size_t cells) {
    return cfd::core::parallel_sum(cells, [&](std::size_t n) {
        double rho = 0.0;
        for (int d = 0; d < 9; ++d) rho += f(static_cast<std::size_t>(d), n);
        return rho;
    });
}

template<class Field>
double kinetic_energy_impl(const Field& f, std::size_t cells) {
    return cfd::core::parallel_sum(cells, [&](std::size_t n) {
        std::array<float, 9> local{};
        for (int d = 0; d < 9; ++d) local[static_cast<std::size_t>(d)] = f(static_cast<std::size_t>(d), n);
        const auto state = macroscopic(local);
        const double ux = state[1];
        const double uy = state[2];
        return 0.5 * static_cast<double>(state[0]) * (ux * ux + uy * uy);
    });
}

constexpr double moment_matrix[9][9] = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1},
    {-4, -1, -1, -1, -1, 2, 2, 2, 2},
    {4, -2, -2, -2, -2, 1, 1, 1, 1},
    {0, 1, 0, -1, 0, 1, -1, -1, 1},
    {0, -2, 0, 2, 0, 1, -1, -1, 1},
    {0, 0, 1, 0, -1, 1, 1, -1, -1},
    {0, 0, -2, 0, 2, 1, 1, -1, -1},
    {0, 1, -1, 1, -1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 1, -1, 1, -1}
};

std::array<std::array<double, 9>, 9> inverse_moment_matrix() {
    double a[9][18]{};
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) a[i][j] = moment_matrix[i][j];
        a[i][9 + i] = 1.0;
    }
    for (int col = 0; col < 9; ++col) {
        int pivot = col;
        for (int row = col + 1; row < 9; ++row) {
            if (std::abs(a[row][col]) > std::abs(a[pivot][col])) pivot = row;
        }
        if (std::abs(a[pivot][col]) < 1.0e-14) throw std::runtime_error("singular MRT transform");
        if (pivot != col) {
            for (int j = 0; j < 18; ++j) std::swap(a[pivot][j], a[col][j]);
        }
        const double diagonal = a[col][col];
        for (int j = 0; j < 18; ++j) a[col][j] /= diagonal;
        for (int row = 0; row < 9; ++row) {
            if (row == col) continue;
            const double factor = a[row][col];
            for (int j = 0; j < 18; ++j) a[row][j] -= factor * a[col][j];
        }
    }
    std::array<std::array<double, 9>, 9> inverse{};
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) inverse[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = a[i][9 + j];
    }
    return inverse;
}

std::array<std::array<double, 9>, 9> inverse_raw_moment_matrix() {
    constexpr int powers[9][2]{{0,0},{1,0},{0,1},{2,0},{0,2},{1,1},{2,1},{1,2},{2,2}};
    double a[9][18]{};
    for(int row=0;row<9;++row){
        for(int d=0;d<9;++d){
            const double cx=static_cast<double>(D2Q9Descriptor::cx(d));
            const double cy=static_cast<double>(D2Q9Descriptor::cy(d));
            const auto ipow=[](double v,int p){return p==0?1.0:(p==1?v:v*v);};
            a[row][d]=ipow(cx,powers[row][0])*ipow(cy,powers[row][1]);
        }
        a[row][9+row]=1.0;
    }
    for(int col=0;col<9;++col){
        int pivot=col;for(int row=col+1;row<9;++row)if(std::abs(a[row][col])>std::abs(a[pivot][col]))pivot=row;
        if(std::abs(a[pivot][col])<1e-14)throw std::runtime_error("singular D2Q9 raw moment transform");
        if(pivot!=col)for(int j=0;j<18;++j)std::swap(a[pivot][j],a[col][j]);
        const double diag=a[col][col];for(int j=0;j<18;++j)a[col][j]/=diag;
        for(int row=0;row<9;++row)if(row!=col){const double factor=a[row][col];for(int j=0;j<18;++j)a[row][j]-=factor*a[col][j];}
    }
    std::array<std::array<double,9>,9> inv{};for(int d=0;d<9;++d)for(int row=0;row<9;++row)inv[static_cast<std::size_t>(d)][static_cast<std::size_t>(row)]=a[d][9+row];return inv;
}

} // namespace

RegularizedD2Q9Solver::RegularizedD2Q9Solver(AdvancedD2Q9Config config)
    : c_(config), grid_(config.nx, config.ny), f_(grid_.cells()), next_(grid_.cells()) {
    if (!(c_.tau > 0.5F) || !std::isfinite(c_.tau)) throw std::invalid_argument("invalid regularized tau");
    initialize_uniform();
}

void RegularizedD2Q9Solver::initialize_uniform(float rho, float ux, float uy) {
    initialize_populations(grid_, f_, [=](std::size_t, std::size_t) {
        return std::array<float, 3>{rho, ux, uy};
    });
}

void RegularizedD2Q9Solver::initialize_taylor_green(float amplitude) {
    const float two_pi = 2.0F * std::numbers::pi_v<float>;
    initialize_populations(grid_, f_, [=, this](std::size_t x, std::size_t y) {
        const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(c_.nx);
        const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(c_.ny);
        return std::array<float, 3>{1.0F,
            amplitude * std::sin(two_pi * xf) * std::cos(two_pi * yf),
            -amplitude * std::cos(two_pi * xf) * std::sin(two_pi * yf)};
    });
}

void RegularizedD2Q9Solver::step_once() {
    const float omega = 1.0F / c_.tau;
    cfd::core::parallel_for(grid_.cells(), [&](std::size_t n) {
        const std::size_t x = n % c_.nx;
        const std::size_t y = n / c_.nx;
        std::array<float, 9> fin{};
        std::array<float, 9> feq{};
        for (int d = 0; d < 9; ++d) {
            fin[static_cast<std::size_t>(d)] = f_(static_cast<std::size_t>(d), periodic_source(grid_, x, y, d));
        }
        const auto state = macroscopic(fin);
        double pxx = 0.0;
        double pxy = 0.0;
        double pyy = 0.0;
        for (int d = 0; d < 9; ++d) {
            const auto i = static_cast<std::size_t>(d);
            feq[i] = D2Q9Descriptor::equilibrium(d, state[0], state[1], state[2]);
            const double neq = static_cast<double>(fin[i] - feq[i]);
            const double cx = static_cast<double>(D2Q9Descriptor::cx(d));
            const double cy = static_cast<double>(D2Q9Descriptor::cy(d));
            pxx += neq * cx * cx;
            pxy += neq * cx * cy;
            pyy += neq * cy * cy;
        }
        for (int d = 0; d < 9; ++d) {
            const auto i = static_cast<std::size_t>(d);
            const double cx = static_cast<double>(D2Q9Descriptor::cx(d));
            const double cy = static_cast<double>(D2Q9Descriptor::cy(d));
            const double regularized = 4.5 * static_cast<double>(D2Q9Descriptor::weight(d)) *
                ((cx * cx - 1.0 / 3.0) * pxx + 2.0 * cx * cy * pxy +
                 (cy * cy - 1.0 / 3.0) * pyy);
            next_(i, n) = static_cast<float>(static_cast<double>(feq[i]) +
                (1.0 - static_cast<double>(omega)) * regularized);
        }
    });
    std::swap(f_, next_);
}

void RegularizedD2Q9Solver::step(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) step_once();
}

double RegularizedD2Q9Solver::mass() const { return mass_impl(f_, grid_.cells()); }
double RegularizedD2Q9Solver::kinetic_energy() const { return kinetic_energy_impl(f_, grid_.cells()); }


CumulantD2Q9Solver::CumulantD2Q9Solver(AdvancedD2Q9Config config)
    : c_(config), grid_(config.nx, config.ny), f_(grid_.cells()), next_(grid_.cells()),
      raw_inv_(inverse_raw_moment_matrix()) {
    if (!(c_.tau > 0.5F) || !std::isfinite(c_.tau)) throw std::invalid_argument("invalid cumulant tau");
    initialize_uniform();
}
void CumulantD2Q9Solver::initialize_uniform(float rho,float ux,float uy){initialize_populations(grid_,f_,[=](std::size_t,std::size_t){return std::array<float,3>{rho,ux,uy};});}
void CumulantD2Q9Solver::initialize_taylor_green(float amplitude){const float two_pi=2.0F*std::numbers::pi_v<float>;initialize_populations(grid_,f_,[=,this](std::size_t x,std::size_t y){const float xf=(static_cast<float>(x)+.5F)/static_cast<float>(c_.nx),yf=(static_cast<float>(y)+.5F)/static_cast<float>(c_.ny);return std::array<float,3>{1.0F,amplitude*std::sin(two_pi*xf)*std::cos(two_pi*yf),-amplitude*std::cos(two_pi*xf)*std::sin(two_pi*yf)};});}
void CumulantD2Q9Solver::step_once(){
    const double omega=1.0/static_cast<double>(c_.tau);constexpr double cs2=1.0/3.0;
    cfd::core::parallel_for(grid_.cells(),[&](std::size_t n){
        const std::size_t x=n%c_.nx,y=n/c_.nx;std::array<float,9>fin{};for(int d=0;d<9;++d)fin[static_cast<std::size_t>(d)]=f_(static_cast<std::size_t>(d),periodic_source(grid_,x,y,d));
        const auto state=macroscopic(fin);const double rho=state[0],ux=state[1],uy=state[2];if(!(rho>0.0)){for(int d=0;d<9;++d)next_(static_cast<std::size_t>(d),n)=fin[static_cast<std::size_t>(d)];return;}
        double mu20=0,mu02=0,mu11=0,mu21=0,mu12=0,mu22=0;
        for(int d=0;d<9;++d){const double dx=static_cast<double>(D2Q9Descriptor::cx(d))-ux,dy=static_cast<double>(D2Q9Descriptor::cy(d))-uy,v=fin[static_cast<std::size_t>(d)];const double dx2=dx*dx,dy2=dy*dy;mu20+=v*dx2;mu02+=v*dy2;mu11+=v*dx*dy;mu21+=v*dx2*dy;mu12+=v*dx*dy2;mu22+=v*dx2*dy2;}
        double k20=mu20/rho,k02=mu02/rho,k11=mu11/rho,k21=mu21/rho,k12=mu12/rho,k22=mu22/rho-k20*k02-2.0*k11*k11;
        k20-=omega*(k20-cs2);k02-=omega*(k02-cs2);k11-=omega*k11;
        // Higher cumulants are non-hydrodynamic on D2Q9; relax them fully to
        // the factorized isothermal equilibrium attractor.
        k21=0.0;k12=0.0;k22=0.0;
        mu20=rho*k20;mu02=rho*k02;mu11=rho*k11;mu21=rho*k21;mu12=rho*k12;mu22=rho*(k22+k20*k02+2.0*k11*k11);
        std::array<double,9> raw{};raw[0]=rho;raw[1]=rho*ux;raw[2]=rho*uy;raw[3]=mu20+rho*ux*ux;raw[4]=mu02+rho*uy*uy;raw[5]=mu11+rho*ux*uy;raw[6]=mu21+2.0*ux*mu11+uy*mu20+rho*ux*ux*uy;raw[7]=mu12+2.0*uy*mu11+ux*mu02+rho*ux*uy*uy;raw[8]=mu22+2.0*uy*mu21+uy*uy*mu20+2.0*ux*mu12+4.0*ux*uy*mu11+ux*ux*mu02+rho*ux*ux*uy*uy;
        for(int d=0;d<9;++d){double value=0;for(int row=0;row<9;++row)value+=raw_inv_[static_cast<std::size_t>(d)][static_cast<std::size_t>(row)]*raw[static_cast<std::size_t>(row)];next_(static_cast<std::size_t>(d),n)=static_cast<float>(value);}
    });std::swap(f_,next_);
}
void CumulantD2Q9Solver::step(std::size_t count){for(std::size_t i=0;i<count;++i)step_once();}
double CumulantD2Q9Solver::mass()const{return mass_impl(f_,grid_.cells());}
double CumulantD2Q9Solver::kinetic_energy()const{return kinetic_energy_impl(f_,grid_.cells());}

MrtD2Q9Solver::MrtD2Q9Solver(AdvancedD2Q9Config config)
    : c_(config), grid_(config.nx, config.ny), f_(grid_.cells()), next_(grid_.cells()),
      inv_(inverse_moment_matrix()) {
    if (!(c_.tau > 0.5F) || !std::isfinite(c_.tau)) throw std::invalid_argument("invalid MRT tau");
    initialize_uniform();
}

void MrtD2Q9Solver::initialize_uniform(float rho, float ux, float uy) {
    initialize_populations(grid_, f_, [=](std::size_t, std::size_t) {
        return std::array<float, 3>{rho, ux, uy};
    });
}

void MrtD2Q9Solver::initialize_taylor_green(float amplitude) {
    const float two_pi = 2.0F * std::numbers::pi_v<float>;
    initialize_populations(grid_, f_, [=, this](std::size_t x, std::size_t y) {
        const float xf = (static_cast<float>(x) + 0.5F) / static_cast<float>(c_.nx);
        const float yf = (static_cast<float>(y) + 0.5F) / static_cast<float>(c_.ny);
        return std::array<float, 3>{1.0F,
            amplitude * std::sin(two_pi * xf) * std::cos(two_pi * yf),
            -amplitude * std::cos(two_pi * xf) * std::sin(two_pi * yf)};
    });
}

void MrtD2Q9Solver::step_once() {
    const double omega = 1.0 / static_cast<double>(c_.tau);
    const std::array<double, 9> relaxation{0.0, 1.1, 1.0, 0.0, 1.2, 0.0, 1.2, omega, omega};
    cfd::core::parallel_for(grid_.cells(), [&](std::size_t n) {
        const std::size_t x = n % c_.nx;
        const std::size_t y = n / c_.nx;
        std::array<float, 9> fin{};
        std::array<float, 9> feq{};
        for (int d = 0; d < 9; ++d) {
            fin[static_cast<std::size_t>(d)] = f_(static_cast<std::size_t>(d), periodic_source(grid_, x, y, d));
        }
        const auto state = macroscopic(fin);
        for (int d = 0; d < 9; ++d) {
            feq[static_cast<std::size_t>(d)] = D2Q9Descriptor::equilibrium(d, state[0], state[1], state[2]);
        }
        std::array<double, 9> moments{};
        std::array<double, 9> equilibrium_moments{};
        for (int row = 0; row < 9; ++row) {
            for (int d = 0; d < 9; ++d) {
                moments[static_cast<std::size_t>(row)] += moment_matrix[row][d] * static_cast<double>(fin[static_cast<std::size_t>(d)]);
                equilibrium_moments[static_cast<std::size_t>(row)] += moment_matrix[row][d] * static_cast<double>(feq[static_cast<std::size_t>(d)]);
            }
            moments[static_cast<std::size_t>(row)] -= relaxation[static_cast<std::size_t>(row)] *
                (moments[static_cast<std::size_t>(row)] - equilibrium_moments[static_cast<std::size_t>(row)]);
        }
        for (int d = 0; d < 9; ++d) {
            double value = 0.0;
            for (int row = 0; row < 9; ++row) {
                value += inv_[static_cast<std::size_t>(d)][static_cast<std::size_t>(row)] * moments[static_cast<std::size_t>(row)];
            }
            next_(static_cast<std::size_t>(d), n) = static_cast<float>(value);
        }
    });
    std::swap(f_, next_);
}

void MrtD2Q9Solver::step(std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) step_once();
}

double MrtD2Q9Solver::mass() const { return mass_impl(f_, grid_.cells()); }
double MrtD2Q9Solver::kinetic_energy() const { return kinetic_energy_impl(f_, grid_.cells()); }

} // namespace cfd::lbm

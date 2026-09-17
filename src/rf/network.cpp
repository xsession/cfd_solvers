#include "cfd/rf/network.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cfd::rf {
namespace {

constexpr double c0 = 299792458.0;

Matrix2C inverse(const Matrix2C& m) {
    const Complex det = m.a11 * m.a22 - m.a12 * m.a21;
    if (std::abs(det) < 1.0e-30) throw std::runtime_error("singular 2x2 RF matrix");
    return {m.a22 / det, -m.a12 / det, -m.a21 / det, m.a11 / det};
}

Matrix2C multiply(const Matrix2C& a, const Matrix2C& b) noexcept {
    return {
        a.a11*b.a11 + a.a12*b.a21,
        a.a11*b.a12 + a.a12*b.a22,
        a.a21*b.a11 + a.a22*b.a21,
        a.a21*b.a12 + a.a22*b.a22
    };
}

Complex touchstone_value(double first, double second, std::string format) {
    std::transform(format.begin(), format.end(), format.begin(), [](unsigned char c){return static_cast<char>(std::toupper(c));});
    if (format == "RI") return {first, second};
    const double angle = second * std::numbers::pi / 180.0;
    if (format == "MA") return std::polar(first, angle);
    if (format == "DB") return std::polar(std::pow(10.0, first / 20.0), angle);
    throw std::invalid_argument("unsupported Touchstone complex format");
}

double unit_scale(std::string unit) {
    std::transform(unit.begin(), unit.end(), unit.begin(), [](unsigned char c){return static_cast<char>(std::toupper(c));});
    if (unit == "HZ") return 1.0;
    if (unit == "KHZ") return 1.0e3;
    if (unit == "MHZ") return 1.0e6;
    if (unit == "GHZ") return 1.0e9;
    throw std::invalid_argument("unsupported Touchstone frequency unit");
}

} // namespace

Matrix2C z_to_s(const Matrix2C& z, double z0) {
    if (!(z0 > 0.0) || !std::isfinite(z0)) throw std::invalid_argument("RF reference impedance must be positive");
    const Matrix2C numerator{z.a11-z0,z.a12,z.a21,z.a22-z0};
    const Matrix2C denominator{z.a11+z0,z.a12,z.a21,z.a22+z0};
    return multiply(numerator, inverse(denominator));
}

Matrix2C s_to_z(const Matrix2C& s, double z0) {
    if (!(z0 > 0.0) || !std::isfinite(z0)) throw std::invalid_argument("RF reference impedance must be positive");
    const Matrix2C plus{1.0+s.a11,s.a12,s.a21,1.0+s.a22};
    const Matrix2C minus{1.0-s.a11,-s.a12,-s.a21,1.0-s.a22};
    const auto result = multiply(plus, inverse(minus));
    return {z0*result.a11,z0*result.a12,z0*result.a21,z0*result.a22};
}

Matrix2C s_to_abcd(const Matrix2C& s, double z0) {
    if (!(z0 > 0.0) || !std::isfinite(z0) || std::abs(s.a21) < 1.0e-30)
        throw std::invalid_argument("invalid S matrix for ABCD conversion");
    const Complex d = 2.0*s.a21;
    return {
        ((1.0+s.a11)*(1.0-s.a22)+s.a12*s.a21)/d,
        z0*((1.0+s.a11)*(1.0+s.a22)-s.a12*s.a21)/d,
        ((1.0-s.a11)*(1.0-s.a22)-s.a12*s.a21)/(z0*d),
        ((1.0-s.a11)*(1.0+s.a22)+s.a12*s.a21)/d
    };
}

Matrix2C abcd_to_s(const Matrix2C& m, double z0) {
    if (!(z0 > 0.0) || !std::isfinite(z0)) throw std::invalid_argument("RF reference impedance must be positive");
    const Complex den = m.a11 + m.a12/z0 + m.a21*z0 + m.a22;
    if (std::abs(den) < 1.0e-30) throw std::runtime_error("singular ABCD to S conversion");
    const Complex determinant = m.a11*m.a22 - m.a12*m.a21;
    return {
        (m.a11 + m.a12/z0 - m.a21*z0 - m.a22)/den,
        2.0*determinant/den,
        2.0/den,
        (-m.a11 + m.a12/z0 - m.a21*z0 + m.a22)/den
    };
}

Matrix2C cascade_abcd(const Matrix2C& lhs, const Matrix2C& rhs) noexcept {
    return multiply(lhs, rhs);
}

Matrix2C cascade_s(const Matrix2C& lhs, const Matrix2C& rhs, double z0) {
    return abcd_to_s(cascade_abcd(s_to_abcd(lhs,z0), s_to_abcd(rhs,z0)), z0);
}

Matrix2C TransmissionLine::abcd() const {
    if (!(characteristic_impedance > 0.0) || !std::isfinite(characteristic_impedance)
        || !(length >= 0.0) || !std::isfinite(length)
        || !std::isfinite(propagation_constant.real()) || !std::isfinite(propagation_constant.imag()))
        throw std::invalid_argument("invalid transmission-line parameters");
    const Complex gl = propagation_constant * length;
    const Complex ch = std::cosh(gl);
    const Complex sh = std::sinh(gl);
    return {ch, characteristic_impedance*sh, sh/characteristic_impedance, ch};
}

Matrix2C TransmissionLine::s_parameters(double z0) const { return abcd_to_s(abcd(), z0); }

TransmissionLine TemTransmissionLine::at_frequency(double frequency) const {
    if (!(characteristic_impedance>0.0) || !(phase_velocity>0.0) || !(attenuation_nepers_per_m>=0.0)
        || !(length>=0.0) || !(frequency>0.0) || !std::isfinite(characteristic_impedance)
        || !std::isfinite(phase_velocity) || !std::isfinite(attenuation_nepers_per_m)
        || !std::isfinite(length) || !std::isfinite(frequency))
        throw std::invalid_argument("invalid TEM transmission-line controls");
    return {characteristic_impedance,
            {attenuation_nepers_per_m,2.0*std::numbers::pi*frequency/phase_velocity},
            length};
}

Matrix2C TemTransmissionLine::s_parameters(double frequency,double z0) const {
    return at_frequency(frequency).s_parameters(z0);
}

MicrostripQuasiStatic microstrip_quasi_static(double width, double height, double er, double frequency) {
    if (!(width > 0.0) || !(height > 0.0) || !(er > 1.0) || !(frequency > 0.0))
        throw std::invalid_argument("invalid microstrip dimensions/material/frequency");
    const double u = width / height;
    const double correction = u < 1.0 ? 0.04*std::pow(1.0-u,2.0) : 0.0;
    const double ee = 0.5*(er+1.0) + 0.5*(er-1.0)*(1.0/std::sqrt(1.0+12.0/u)+correction);
    double z0 = 0.0;
    if (u <= 1.0) {
        z0 = (60.0/std::sqrt(ee))*std::log(8.0/u + 0.25*u);
    } else {
        z0 = (120.0*std::numbers::pi)/(std::sqrt(ee)*(u+1.393+0.667*std::log(u+1.444)));
    }
    const double vp = c0/std::sqrt(ee);
    return {z0,ee,vp,vp/frequency};
}

CoaxQuasiStatic coax_quasi_static(double inner_radius,double outer_radius,double er,double frequency) {
    if (!(inner_radius>0.0) || !(outer_radius>inner_radius) || !(er>0.0) || !(frequency>0.0)
        || !std::isfinite(inner_radius) || !std::isfinite(outer_radius) || !std::isfinite(er) || !std::isfinite(frequency))
        throw std::invalid_argument("invalid coax dimensions/material/frequency");
    constexpr double eta0 = 376.730313668;
    const double velocity=c0/std::sqrt(er);
    const double impedance=eta0/(2.0*std::numbers::pi*std::sqrt(er))*std::log(outer_radius/inner_radius);
    const double beta=2.0*std::numbers::pi*frequency/velocity;
    return {impedance,velocity,beta,velocity/frequency};
}

RectangularWaveguideTE10 rectangular_waveguide_te10(double a,double b,double er,double mur,double frequency) {
    if (!(a>0.0) || !(b>0.0) || !(er>0.0) || !(mur>0.0) || !(frequency>0.0)
        || !std::isfinite(a) || !std::isfinite(b) || !std::isfinite(er) || !std::isfinite(mur) || !std::isfinite(frequency))
        throw std::invalid_argument("invalid rectangular waveguide dimensions/material/frequency");
    (void)b; // TE10 cutoff depends only on the broad wall.
    constexpr double eta0 = 376.730313668;
    const double medium_factor=std::sqrt(er*mur);
    const double cutoff=c0/(2.0*a*medium_factor);
    if (!(frequency>cutoff)) throw std::domain_error("TE10 frequency must be above cutoff");
    const double ratio=cutoff/frequency;
    const double root=std::sqrt(1.0-ratio*ratio);
    const double beta=2.0*std::numbers::pi*frequency*medium_factor/c0*root;
    const double wavelength=2.0*std::numbers::pi/beta;
    const double impedance=eta0*std::sqrt(mur/er)/root;
    return {cutoff,beta,wavelength,impedance};
}

std::vector<Touchstone2PortPoint> read_touchstone_s2p(std::istream& input, double* z0_out) {
    std::string unit="GHZ", parameter="S", format="MA";
    double z0=50.0;
    std::vector<double> pending;
    std::vector<Touchstone2PortPoint> points;
    std::string line;
    while (std::getline(input,line)) {
        const auto bang=line.find('!'); if(bang!=std::string::npos) line.erase(bang);
        if(line.empty()) continue;
        if(line.front()=='#') {
            std::istringstream ss(line.substr(1)); ss>>unit>>parameter>>format;
            std::string token; while(ss>>token){
                std::string upper=token;std::transform(upper.begin(),upper.end(),upper.begin(),[](unsigned char c){return static_cast<char>(std::toupper(c));});
                if(upper=="R") ss>>z0;
            }
            std::transform(parameter.begin(),parameter.end(),parameter.begin(),[](unsigned char c){return static_cast<char>(std::toupper(c));});
            if(parameter!="S") throw std::invalid_argument("Touchstone reader expects S-parameters");
            continue;
        }
        std::istringstream ss(line); double value=0.0; while(ss>>value) pending.push_back(value);
        while(pending.size()>=9U){
            const double scale=unit_scale(unit);
            Matrix2C s;
            s.a11=touchstone_value(pending[1],pending[2],format);
            s.a21=touchstone_value(pending[3],pending[4],format);
            s.a12=touchstone_value(pending[5],pending[6],format);
            s.a22=touchstone_value(pending[7],pending[8],format);
            points.push_back({pending[0]*scale,s});
            pending.erase(pending.begin(),pending.begin()+9);
        }
    }
    if(!pending.empty()) throw std::invalid_argument("truncated Touchstone S2P data");
    if(z0_out) *z0_out=z0;
    return points;
}

void write_touchstone_s2p(std::ostream& output, const std::vector<Touchstone2PortPoint>& points, double z0) {
    if (!(z0 > 0.0) || !std::isfinite(z0)) throw std::invalid_argument("invalid Touchstone reference impedance");
    output << "# HZ S RI R " << z0 << '\n';
    for (const auto& p:points) {
        output << p.frequency_hz << ' '
               << p.s.a11.real() << ' ' << p.s.a11.imag() << ' '
               << p.s.a21.real() << ' ' << p.s.a21.imag() << ' '
               << p.s.a12.real() << ' ' << p.s.a12.imag() << ' '
               << p.s.a22.real() << ' ' << p.s.a22.imag() << '\n';
    }
}


PlanarLineQuasiStatic stripline_quasi_static(double width,double spacing,double er,double frequency) {
    if(!(width>0.0)||!(spacing>0.0)||!(er>0.0)||!(frequency>0.0)
       ||!std::isfinite(width)||!std::isfinite(spacing)||!std::isfinite(er)||!std::isfinite(frequency))
        throw std::invalid_argument("invalid stripline parameters");
    // Wheeler-style thin-conductor symmetric stripline approximation.
    const double u=width/spacing;
    const double z0=(30.0*std::numbers::pi/std::sqrt(er))/(u+0.441);
    const double velocity=299792458.0/std::sqrt(er);
    return {z0,er,velocity,velocity/frequency};
}

namespace {
double complete_elliptic_k(double k) {
    if(!(k>=0.0&&k<1.0)) throw std::invalid_argument("elliptic modulus must satisfy 0 <= k < 1");
    double a=1.0,b=std::sqrt(1.0-k*k);
    for(int i=0;i<32 && std::abs(a-b)>1.0e-15*a;++i){const double next=0.5*(a+b);b=std::sqrt(a*b);a=next;}
    return std::numbers::pi/(2.0*a);
}
}

PlanarLineQuasiStatic coplanar_waveguide_quasi_static(double width,double gap,double er,double frequency) {
    if(!(width>0.0)||!(gap>0.0)||!(er>0.0)||!(frequency>0.0)
       ||!std::isfinite(width)||!std::isfinite(gap)||!std::isfinite(er)||!std::isfinite(frequency))
        throw std::invalid_argument("invalid coplanar waveguide parameters");
    const double k=width/(width+2.0*gap);
    const double kp=std::sqrt(1.0-k*k);
    const double eeff=0.5*(er+1.0);
    const double z0=30.0*std::numbers::pi/std::sqrt(eeff)*complete_elliptic_k(kp)/complete_elliptic_k(k);
    const double velocity=299792458.0/std::sqrt(eeff);
    return {z0,eeff,velocity,velocity/frequency};
}

} // namespace cfd::rf

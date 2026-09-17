#include "cfd/rf/nport.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cfd::rf {
namespace {

std::string upper(std::string text) {
    for (char& c : text) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return text;
}

double frequency_scale(const std::string& input) {
    const std::string unit = upper(input);
    if (unit == "HZ") return 1.0;
    if (unit == "KHZ") return 1.0e3;
    if (unit == "MHZ") return 1.0e6;
    if (unit == "GHZ") return 1.0e9;
    throw std::invalid_argument("unsupported Touchstone frequency unit");
}

Complex decode_pair(double first, double second, const std::string& input_format) {
    const std::string format = upper(input_format);
    if (format == "RI") return {first, second};
    const double radians = second * std::numbers::pi / 180.0;
    if (format == "MA") return std::polar(first, radians);
    if (format == "DB") return std::polar(std::pow(10.0, first / 20.0), radians);
    throw std::invalid_argument("unsupported Touchstone data format");
}

void validate_reference(std::span<const double> reference, std::size_t n) {
    if (reference.size() != n) throw std::invalid_argument("N-port reference-impedance size mismatch");
    for (double z0 : reference) {
        if (!(z0 > 0.0) || !std::isfinite(z0)) throw std::invalid_argument("N-port reference impedance must be finite and positive");
    }
}

void validate_power_reference(std::span<const Complex> reference, std::size_t n) {
    if (reference.size() != n) throw std::invalid_argument("N-port complex reference-impedance size mismatch");
    for (Complex z0 : reference) {
        if (!(z0.real() > 0.0) || !std::isfinite(z0.real()) || !std::isfinite(z0.imag()))
            throw std::invalid_argument("power-wave reference impedance must be finite with positive real part");
    }
}

ComplexMatrix diagonal_real_sqrt(std::span<const Complex> reference, bool reciprocal) {
    ComplexMatrix result(reference.size());
    for (std::size_t i = 0; i < reference.size(); ++i) {
        const double root = std::sqrt(reference[i].real());
        result(i, i) = reciprocal ? 1.0 / root : root;
    }
    return result;
}

ComplexMatrix diagonal_complex(std::span<const Complex> reference, bool conjugate) {
    ComplexMatrix result(reference.size());
    for (std::size_t i = 0; i < reference.size(); ++i)
        result(i, i) = conjugate ? std::conj(reference[i]) : reference[i];
    return result;
}

ComplexMatrix diagonal_sqrt(std::span<const double> reference, bool reciprocal) {
    ComplexMatrix result(reference.size());
    for (std::size_t i = 0; i < reference.size(); ++i) {
        const double root = std::sqrt(reference[i]);
        result(i, i) = reciprocal ? 1.0 / root : root;
    }
    return result;
}

ComplexMatrix add_identity(const ComplexMatrix& matrix, double sign) {
    ComplexMatrix result = matrix;
    for (std::size_t i = 0; i < matrix.size(); ++i) result(i, i) += sign;
    return result;
}


std::vector<double> parse_numbers(std::string line) {
    const auto comment = line.find('!');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream stream(line);
    std::vector<double> values;
    double value = 0.0;
    while (stream >> value) values.push_back(value);
    return values;
}

} // namespace

ComplexMatrix::ComplexMatrix(std::size_t size, Complex value)
    : size_(size), data_(size * size, value) {}

Complex& ComplexMatrix::operator()(std::size_t row, std::size_t column) {
    if (row >= size_ || column >= size_) throw std::out_of_range("RF matrix index out of range");
    return data_[row * size_ + column];
}

const Complex& ComplexMatrix::operator()(std::size_t row, std::size_t column) const {
    if (row >= size_ || column >= size_) throw std::out_of_range("RF matrix index out of range");
    return data_[row * size_ + column];
}

ComplexMatrix ComplexMatrix::identity(std::size_t size) {
    ComplexMatrix result(size);
    for (std::size_t i = 0; i < size; ++i) result(i, i) = 1.0;
    return result;
}

ComplexMatrix multiply(const ComplexMatrix& lhs, const ComplexMatrix& rhs) {
    if (lhs.size() != rhs.size()) throw std::invalid_argument("RF matrix multiply size mismatch");
    const std::size_t n = lhs.size();
    ComplexMatrix result(n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t k = 0; k < n; ++k) {
            const Complex a = lhs(row, k);
            for (std::size_t column = 0; column < n; ++column) result(row, column) += a * rhs(k, column);
        }
    }
    return result;
}

ComplexMatrix inverse(const ComplexMatrix& matrix) {
    const std::size_t n = matrix.size();
    if (n == 0U) throw std::invalid_argument("cannot invert empty RF matrix");
    std::vector<Complex> augmented(n * 2U * n, Complex{});
    const std::size_t stride = 2U * n;
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t column = 0; column < n; ++column) augmented[row * stride + column] = matrix(row, column);
        augmented[row * stride + n + row] = 1.0;
    }
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        double best = std::abs(augmented[column * stride + column]);
        for (std::size_t row = column + 1U; row < n; ++row) {
            const double value = std::abs(augmented[row * stride + column]);
            if (value > best) { best = value; pivot = row; }
        }
        if (best < 1.0e-30) throw std::runtime_error("singular RF N-port matrix");
        if (pivot != column) {
            for (std::size_t j = 0; j < stride; ++j) std::swap(augmented[pivot * stride + j], augmented[column * stride + j]);
        }
        const Complex diagonal = augmented[column * stride + column];
        for (std::size_t j = 0; j < stride; ++j) augmented[column * stride + j] /= diagonal;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == column) continue;
            const Complex factor = augmented[row * stride + column];
            if (std::abs(factor) == 0.0) continue;
            for (std::size_t j = 0; j < stride; ++j) augmented[row * stride + j] -= factor * augmented[column * stride + j];
        }
    }
    ComplexMatrix result(n);
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column)
            result(row, column) = augmented[row * stride + n + column];
    return result;
}

ComplexMatrix z_to_s(const ComplexMatrix& z, std::span<const double> reference) {
    const std::size_t n = z.size();
    validate_reference(reference, n);
    const ComplexMatrix root = diagonal_sqrt(reference, false);
    const ComplexMatrix inverse_root = diagonal_sqrt(reference, true);
    // Normalize Z to z_n = Z0^-1/2 Z Z0^-1/2, then S=(z_n-I)(z_n+I)^-1.
    const ComplexMatrix normalized = multiply(multiply(inverse_root, z), inverse_root);
    ComplexMatrix minus = normalized;
    ComplexMatrix plus = normalized;
    for (std::size_t i = 0; i < n; ++i) { minus(i, i) -= 1.0; plus(i, i) += 1.0; }
    return multiply(minus, inverse(plus));
}

ComplexMatrix s_to_z(const ComplexMatrix& s, std::span<const double> reference) {
    const std::size_t n = s.size();
    validate_reference(reference, n);
    const ComplexMatrix root = diagonal_sqrt(reference, false);
    const ComplexMatrix plus = add_identity(s, 1.0);
    ComplexMatrix minus = ComplexMatrix::identity(n);
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column)
            minus(row, column) -= s(row, column);
    const ComplexMatrix normalized = multiply(plus, inverse(minus));
    return multiply(multiply(root, normalized), root);
}

ComplexMatrix renormalize_s(const ComplexMatrix& s, std::span<const double> old_reference,
                            std::span<const double> new_reference) {
    return z_to_s(s_to_z(s, old_reference), new_reference);
}

ComplexMatrix z_to_s_power_wave(const ComplexMatrix& z, std::span<const Complex> reference) {
    const std::size_t n = z.size();
    validate_power_reference(reference, n);
    const ComplexMatrix r_root = diagonal_real_sqrt(reference, false);
    const ComplexMatrix r_inverse_root = diagonal_real_sqrt(reference, true);
    const ComplexMatrix z0 = diagonal_complex(reference, false);
    const ComplexMatrix z0_conjugate = diagonal_complex(reference, true);
    ComplexMatrix numerator = z;
    ComplexMatrix denominator = z;
    for (std::size_t i = 0; i < n; ++i) {
        numerator(i, i) -= z0_conjugate(i, i);
        denominator(i, i) += z0(i, i);
    }
    // b = R^-1/2 (Z-Z0*) (Z+Z0)^-1 R^1/2 a
    return multiply(multiply(multiply(r_inverse_root, numerator), inverse(denominator)), r_root);
}

ComplexMatrix s_to_z_power_wave(const ComplexMatrix& s, std::span<const Complex> reference) {
    const std::size_t n = s.size();
    validate_power_reference(reference, n);
    const ComplexMatrix r_root = diagonal_real_sqrt(reference, false);
    const ComplexMatrix r_inverse_root = diagonal_real_sqrt(reference, true);
    const ComplexMatrix z0 = diagonal_complex(reference, false);
    const ComplexMatrix z0_conjugate = diagonal_complex(reference, true);
    // T = R^1/2 S R^-1/2 and (I-T) Z = T Z0 + Z0*.
    const ComplexMatrix t = multiply(multiply(r_root, s), r_inverse_root);
    ComplexMatrix lhs = ComplexMatrix::identity(n);
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column)
            lhs(row, column) -= t(row, column);
    ComplexMatrix rhs = multiply(t, z0);
    for (std::size_t i = 0; i < n; ++i) rhs(i, i) += z0_conjugate(i, i);
    return multiply(inverse(lhs), rhs);
}

ComplexMatrix renormalize_s_power_wave(const ComplexMatrix& s, std::span<const Complex> old_reference,
                                       std::span<const Complex> new_reference) {
    return z_to_s_power_wave(s_to_z_power_wave(s, old_reference), new_reference);
}

std::vector<NPortPoint> read_touchstone(std::istream& input, std::size_t port_hint) {
    std::string unit = "GHZ";
    std::string format = "MA";
    std::size_t ports = port_hint;
    double scalar_reference = 50.0;
    std::vector<double> explicit_reference;
    std::vector<double> pending;
    std::vector<NPortPoint> points;
    bool in_network_data = false;
    std::string line;
    while (std::getline(input, line)) {
        const auto comment = line.find('!');
        if (comment != std::string::npos) line.erase(comment);
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        line.erase(0, first);
        if (line.front() == '#') {
            std::istringstream options(line.substr(1));
            std::string parameter;
            options >> unit >> parameter >> format;
            if (upper(parameter) != "S") throw std::invalid_argument("Touchstone N-port reader supports S parameters only");
            std::string token;
            while (options >> token) if (upper(token) == "R") options >> scalar_reference;
            continue;
        }
        if (line.front() == '[') {
            const auto close = line.find(']');
            if (close == std::string::npos) throw std::invalid_argument("malformed Touchstone keyword");
            const std::string keyword = upper(line.substr(1, close - 1));
            const std::string rest = line.substr(close + 1);
            if (keyword == "NUMBER OF PORTS") {
                std::istringstream(rest) >> ports;
            } else if (keyword == "REFERENCE") {
                explicit_reference = parse_numbers(rest);
            } else if (keyword == "NETWORK DATA") {
                in_network_data = true;
            } else if (keyword == "END") {
                break;
            }
            continue;
        }
        (void)in_network_data;
        const auto values = parse_numbers(line);
        pending.insert(pending.end(), values.begin(), values.end());
        if (ports == 0U) {
            // v1 SnP needs a filename to infer N. Without one, accept S2P for compatibility.
            ports = 2U;
        }
        const std::size_t count = 1U + 2U * ports * ports;
        while (pending.size() >= count) {
            NPortPoint point;
            point.frequency_hz = pending[0] * frequency_scale(unit);
            point.s = ComplexMatrix(ports);
            point.reference_impedance.assign(ports, scalar_reference);
            if (!explicit_reference.empty()) {
                if (explicit_reference.size() == 1U) point.reference_impedance.assign(ports, explicit_reference.front());
                else {
                    if (explicit_reference.size() != ports) throw std::invalid_argument("Touchstone [Reference] count mismatch");
                    point.reference_impedance = explicit_reference;
                }
            }
            std::size_t cursor = 1U;
            for (std::size_t column = 0; column < ports; ++column) {
                for (std::size_t row = 0; row < ports; ++row) {
                    point.s(row, column) = decode_pair(pending[cursor], pending[cursor + 1U], format);
                    cursor += 2U;
                }
            }
            points.push_back(std::move(point));
            pending.erase(pending.begin(), pending.begin() + static_cast<std::ptrdiff_t>(count));
        }
    }
    if (!pending.empty()) throw std::invalid_argument("truncated Touchstone SnP data");
    return points;
}

void write_touchstone(std::ostream& output, std::span<const NPortPoint> points, bool version2) {
    if (points.empty()) return;
    const std::size_t ports = points.front().s.size();
    if (ports == 0U) throw std::invalid_argument("cannot write empty-port Touchstone network");
    for (const auto& point : points) {
        if (point.s.size() != ports || point.reference_impedance.size() != ports)
            throw std::invalid_argument("Touchstone point topology mismatch");
        validate_reference(point.reference_impedance, ports);
    }
    if (version2) {
        output << "[Version] 2.0\n[Number of Ports] " << ports << "\n";
        output << "[Reference]";
        for (double z0 : points.front().reference_impedance) output << ' ' << z0;
        output << "\n# HZ S RI R 50\n[Network Data]\n";
    } else {
        const double z0 = points.front().reference_impedance.front();
        for (double value : points.front().reference_impedance)
            if (std::abs(value - z0) > 1.0e-12) throw std::invalid_argument("Touchstone v1 writer requires one common reference impedance");
        output << "# HZ S RI R " << z0 << '\n';
    }
    for (const auto& point : points) {
        output << point.frequency_hz;
        for (std::size_t column = 0; column < ports; ++column)
            for (std::size_t row = 0; row < ports; ++row) {
                const Complex value = point.s(row, column);
                output << ' ' << value.real() << ' ' << value.imag();
            }
        output << '\n';
    }
    if (version2) output << "[End]\n";
}

ComplexMatrix shift_reference_planes(const ComplexMatrix& s, std::span<const Complex> gamma,
                                     std::span<const double> distance) {
    const std::size_t n = s.size();
    if (gamma.size() != n || distance.size() != n) throw std::invalid_argument("reference-plane shift size mismatch");
    ComplexMatrix result(n);
    std::vector<Complex> factor(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(distance[i]) || !std::isfinite(gamma[i].real()) || !std::isfinite(gamma[i].imag()))
            throw std::invalid_argument("invalid reference-plane shift");
        factor[i] = std::exp(gamma[i] * distance[i]);
    }
    for (std::size_t row = 0; row < n; ++row)
        for (std::size_t column = 0; column < n; ++column)
            result(row, column) = s(row, column) * factor[row] * factor[column];
    return result;
}

double return_loss_db(Complex reflection) {
    const double magnitude = std::abs(reflection);
    if (!std::isfinite(magnitude)) throw std::invalid_argument("invalid reflection coefficient");
    if (magnitude == 0.0) return std::numeric_limits<double>::infinity();
    return -20.0 * std::log10(magnitude);
}

double vswr(Complex reflection) {
    const double magnitude = std::abs(reflection);
    if (!(magnitude >= 0.0) || !std::isfinite(magnitude)) throw std::invalid_argument("invalid reflection coefficient");
    if (magnitude >= 1.0) return std::numeric_limits<double>::infinity();
    return (1.0 + magnitude) / (1.0 - magnitude);
}

StabilityMetrics stability_metrics(const Matrix2C& s) {
    const Complex delta = s.a11 * s.a22 - s.a12 * s.a21;
    const double denominator = 2.0 * std::abs(s.a12 * s.a21);
    const double k = denominator > 0.0
        ? (1.0 - std::norm(s.a11) - std::norm(s.a22) + std::norm(delta)) / denominator
        : std::numeric_limits<double>::infinity();
    const double mu_den = std::abs(s.a22 - std::conj(s.a11) * delta) + std::abs(s.a12 * s.a21);
    const double mu = mu_den > 0.0 ? (1.0 - std::norm(s.a11)) / mu_den : std::numeric_limits<double>::infinity();
    return {k, mu, delta};
}



StabilityCircle source_stability_circle(const Matrix2C& s) {
    const Complex delta=s.a11*s.a22-s.a12*s.a21;
    const double denominator=std::norm(s.a11)-std::norm(delta);
    if(std::abs(denominator)<1.0e-30)return {{},std::numeric_limits<double>::infinity(),false};
    const Complex center=std::conj(s.a11-delta*std::conj(s.a22))/denominator;
    const double radius=std::abs(s.a12*s.a21)/std::abs(denominator);
    return {center,radius,std::isfinite(radius)&&std::isfinite(center.real())&&std::isfinite(center.imag())};
}

StabilityCircle load_stability_circle(const Matrix2C& s) {
    const Complex delta=s.a11*s.a22-s.a12*s.a21;
    const double denominator=std::norm(s.a22)-std::norm(delta);
    if(std::abs(denominator)<1.0e-30)return {{},std::numeric_limits<double>::infinity(),false};
    const Complex center=std::conj(s.a22-delta*std::conj(s.a11))/denominator;
    const double radius=std::abs(s.a12*s.a21)/std::abs(denominator);
    return {center,radius,std::isfinite(radius)&&std::isfinite(center.real())&&std::isfinite(center.imag())};
}

double transducer_gain(const Matrix2C& s,Complex gamma_s,Complex gamma_l) {
    if(std::abs(gamma_s)>=1.0||std::abs(gamma_l)>=1.0)throw std::invalid_argument("source/load reflection magnitude must be < 1");
    const Complex denominator=(Complex{1.0,0.0}-s.a11*gamma_s)*(Complex{1.0,0.0}-s.a22*gamma_l)
        -s.a12*s.a21*gamma_s*gamma_l;
    const double d2=std::norm(denominator);
    if(d2<1.0e-30)return std::numeric_limits<double>::infinity();
    return (1.0-std::norm(gamma_s))*std::norm(s.a21)*(1.0-std::norm(gamma_l))/d2;
}

double noise_factor(const NoiseParameters& parameters,Complex gamma_s) {
    if(!(parameters.minimum_noise_factor>=1.0)||!(parameters.equivalent_noise_resistance_ohm>=0.0)
       ||!(parameters.reference_impedance_ohm>0.0)||!std::isfinite(parameters.minimum_noise_factor)
       ||!std::isfinite(parameters.equivalent_noise_resistance_ohm)||!std::isfinite(parameters.reference_impedance_ohm)
       ||std::abs(parameters.optimum_source_reflection)>=1.0||std::abs(gamma_s)>=1.0)
        throw std::invalid_argument("invalid RF noise parameters/reflection coefficient");
    const double denominator=(1.0-std::norm(gamma_s))*std::norm(Complex{1.0,0.0}+parameters.optimum_source_reflection);
    const double coefficient=4.0*parameters.equivalent_noise_resistance_ohm/parameters.reference_impedance_ohm;
    return parameters.minimum_noise_factor+coefficient*std::norm(gamma_s-parameters.optimum_source_reflection)/denominator;
}

NoiseCircle noise_circle(const NoiseParameters& parameters,double target) {
    if(!(target>=parameters.minimum_noise_factor)||!std::isfinite(target)
       ||!(parameters.equivalent_noise_resistance_ohm>0.0)||!(parameters.reference_impedance_ohm>0.0)
       ||std::abs(parameters.optimum_source_reflection)>=1.0)
        throw std::invalid_argument("invalid RF noise-circle controls");
    const double coefficient=4.0*parameters.equivalent_noise_resistance_ohm/parameters.reference_impedance_ohm;
    const double n=(target-parameters.minimum_noise_factor)
        *std::norm(Complex{1.0,0.0}+parameters.optimum_source_reflection)/coefficient;
    const double denominator=1.0+n;
    const Complex center=parameters.optimum_source_reflection/denominator;
    const double radius_squared=n*(denominator-std::norm(parameters.optimum_source_reflection))/(denominator*denominator);
    const double radius=std::sqrt(std::max(0.0,radius_squared));
    return {center,radius,target,std::isfinite(radius)&&std::abs(center)+radius<1.0+1.0e-12};
}

std::vector<LoadPullSample> sample_load_pull(const Matrix2C& s,Complex gamma_s,std::size_t radial_steps,
                                              std::size_t angular_steps,double maximum_radius) {
    if(radial_steps==0U||angular_steps==0U||!(maximum_radius>=0.0&&maximum_radius<1.0))
        throw std::invalid_argument("invalid load-pull sampling controls");
    std::vector<LoadPullSample> result;result.reserve(1U+radial_steps*angular_steps);
    result.push_back({{},transducer_gain(s,gamma_s,{})});
    for(std::size_t r=1;r<=radial_steps;++r){
        const double radius=maximum_radius*static_cast<double>(r)/static_cast<double>(radial_steps);
        for(std::size_t a=0;a<angular_steps;++a){
            const double angle=2.0*std::numbers::pi*static_cast<double>(a)/static_cast<double>(angular_steps);
            const Complex gamma=std::polar(radius,angle);
            result.push_back({gamma,transducer_gain(s,gamma_s,gamma)});
        }
    }
    return result;
}

NetworkQuality network_quality(const ComplexMatrix& s,double passivity_tolerance,double reciprocity_tolerance) {
    const std::size_t n=s.size();
    if(n==0U || !(passivity_tolerance>=0.0) || !(reciprocity_tolerance>=0.0))
        throw std::invalid_argument("invalid network-quality controls");
    // Power iteration on S^H S gives sigma_max^2 without a heavyweight SVD dependency.
    std::vector<Complex> x(n,Complex{1.0/std::sqrt(static_cast<double>(n)),0.0}),y(n),z(n);
    double lambda=0.0;
    for(std::size_t iteration=0;iteration<80U;++iteration){
        std::fill(y.begin(),y.end(),Complex{});
        std::fill(z.begin(),z.end(),Complex{});
        for(std::size_t r=0;r<n;++r)for(std::size_t c=0;c<n;++c)y[r]+=s(r,c)*x[c];
        for(std::size_t c=0;c<n;++c)for(std::size_t r=0;r<n;++r)z[c]+=std::conj(s(r,c))*y[r];
        double norm2=0.0;for(const auto& value:z)norm2+=std::norm(value);
        const double norm=std::sqrt(norm2);if(norm<1.0e-30){lambda=0.0;break;}
        for(std::size_t i=0;i<n;++i)x[i]=z[i]/norm;
        Complex rq{};for(std::size_t i=0;i<n;++i)rq+=std::conj(x[i])*z[i];
        const double next=std::max(0.0,rq.real());
        if(std::abs(next-lambda)<1.0e-14*std::max(1.0,next)){lambda=next;break;}
        lambda=next;
    }
    const double sigma=std::sqrt(std::max(0.0,lambda));
    double diff2=0.0,scale2=0.0;
    for(std::size_t r=0;r<n;++r)for(std::size_t c=0;c<n;++c){
        diff2+=std::norm(s(r,c)-s(c,r));scale2+=std::norm(s(r,c));
    }
    const double reciprocity=std::sqrt(diff2)/std::max(std::sqrt(scale2),1.0e-30);
    return {sigma,reciprocity,sigma<=1.0+passivity_tolerance,reciprocity<=reciprocity_tolerance};
}

ComplexMatrix enforce_passivity(const ComplexMatrix& s,double maximum_singular_value) {
    if(s.size()==0U||!(maximum_singular_value>0.0)||!std::isfinite(maximum_singular_value))
        throw std::invalid_argument("invalid passivity-enforcement controls");
    const double sigma=network_quality(s).maximum_singular_value;
    if(!std::isfinite(sigma))throw std::runtime_error("non-finite network singular value");
    if(sigma<=maximum_singular_value)return s;
    ComplexMatrix projected=s;
    const double scale=maximum_singular_value/sigma;
    for(std::size_t row=0;row<s.size();++row)for(std::size_t column=0;column<s.size();++column)
        projected(row,column)*=scale;
    return projected;
}

BroadbandCausality check_broadband_causality(std::span<const NPortPoint> points,
                                              double tolerance,double frequency_tolerance) {
    if(points.size()<3U||!(tolerance>=0.0)||!std::isfinite(tolerance)
       ||!(frequency_tolerance>=0.0)||!std::isfinite(frequency_tolerance))
        throw std::invalid_argument("invalid broadband causality controls");
    const std::size_t ports=points.front().s.size();
    if(ports==0U)throw std::invalid_argument("causality check requires non-empty networks");
    if(std::abs(points.front().frequency_hz)>frequency_tolerance)
        throw std::invalid_argument("causality check requires a DC first sample");
    const double df=points[1].frequency_hz-points[0].frequency_hz;
    if(!(df>0.0)||!std::isfinite(df))throw std::invalid_argument("causality check requires increasing frequencies");
    for(std::size_t k=0;k<points.size();++k){
        if(points[k].s.size()!=ports)throw std::invalid_argument("causality sweep port-count mismatch");
        const double expected=df*static_cast<double>(k);
        const double scale=std::max({1.0,std::abs(expected),std::abs(points[k].frequency_hz)});
        if(std::abs(points[k].frequency_hz-expected)>frequency_tolerance*scale)
            throw std::invalid_argument("causality check requires uniformly spaced frequencies");
    }
    const std::size_t spectrum_size=2U*(points.size()-1U);
    const std::size_t first_negative=points.size();
    double negative_energy=0.0,total_energy=0.0;
    std::vector<Complex> spectrum(spectrum_size),impulse(spectrum_size);
    for(std::size_t row=0;row<ports;++row)for(std::size_t column=0;column<ports;++column){
        std::fill(spectrum.begin(),spectrum.end(),Complex{});
        for(std::size_t k=0;k<points.size();++k)spectrum[k]=points[k].s(row,column);
        // A real-valued impulse response requires real DC/Nyquist bins. Small
        // numerical imaginary residue is discarded at those self-conjugate bins.
        spectrum[0]=Complex{spectrum[0].real(),0.0};
        spectrum[points.size()-1U]=Complex{spectrum[points.size()-1U].real(),0.0};
        for(std::size_t k=1U;k+1U<points.size();++k)spectrum[spectrum_size-k]=std::conj(spectrum[k]);
        for(std::size_t n=0;n<spectrum_size;++n){
            Complex value{};
            for(std::size_t k=0;k<spectrum_size;++k){
                const double phase=2.0*std::numbers::pi*static_cast<double>(k*n)/static_cast<double>(spectrum_size);
                value+=spectrum[k]*std::exp(Complex{0.0,phase});
            }
            impulse[n]=value/static_cast<double>(spectrum_size);
            const double energy=std::norm(impulse[n]);total_energy+=energy;
            if(n>=first_negative)negative_energy+=energy;
        }
    }
    const double ratio=total_energy>0.0?negative_energy/total_energy:0.0;
    return {ratio,ratio<=tolerance,spectrum_size};
}

AdaptiveNetworkSweepResult adaptive_network_sweep(double start,double stop,
                                                         std::span<const double> references,
                                                         NetworkSampler sampler,
                                                         const AdaptiveNetworkSweepConfig& config) {
    if(!(start>0.0)||!(stop>start)||references.empty()||!sampler||config.initial_points<2U
       ||config.maximum_points<config.initial_points||config.maximum_refinements==0U
       ||!(config.relative_tolerance>0.0)||!(config.absolute_tolerance>0.0))
        throw std::invalid_argument("invalid adaptive network sweep controls");
    for(double z0:references)if(!(z0>0.0)||!std::isfinite(z0))throw std::invalid_argument("invalid network reference impedance");
    AdaptiveNetworkSweepResult result;
    const auto sample=[&](double frequency){
        NPortPoint point;point.frequency_hz=frequency;point.s=sampler(frequency);++result.evaluations;
        if(point.s.size()!=references.size())throw std::runtime_error("adaptive network sampler port-count mismatch");
        point.reference_impedance.assign(references.begin(),references.end());return point;
    };
    result.points.reserve(config.maximum_points);
    const double ratio=std::pow(stop/start,1.0/static_cast<double>(config.initial_points-1U));
    for(std::size_t i=0;i<config.initial_points;++i)
        result.points.push_back(sample(i+1U==config.initial_points?stop:start*std::pow(ratio,static_cast<double>(i))));
    for(std::size_t refinement=0;refinement<config.maximum_refinements;++refinement){
        std::vector<NPortPoint> additions;
        for(std::size_t i=0;i+1U<result.points.size();++i){
            if(result.points.size()+additions.size()>=config.maximum_points)break;
            const auto& left=result.points[i];const auto& right=result.points[i+1U];
            const double midpoint=std::sqrt(left.frequency_hz*right.frequency_hz);
            auto actual=sample(midpoint);double normalized_error=0.0;
            for(std::size_t r=0;r<actual.s.size();++r)for(std::size_t c=0;c<actual.s.size();++c){
                const Complex predicted=0.5*(left.s(r,c)+right.s(r,c));
                const double scale=config.absolute_tolerance+config.relative_tolerance*std::max({std::abs(actual.s(r,c)),std::abs(left.s(r,c)),std::abs(right.s(r,c))});
                normalized_error=std::max(normalized_error,std::abs(actual.s(r,c)-predicted)/scale);
            }
            if(normalized_error>1.0)additions.push_back(std::move(actual));
        }
        if(additions.empty()){result.converged=true;break;}
        ++result.refinements;
        result.points.insert(result.points.end(),std::make_move_iterator(additions.begin()),std::make_move_iterator(additions.end()));
        std::sort(result.points.begin(),result.points.end(),[](const auto& a,const auto& b){return a.frequency_hz<b.frequency_hz;});
        if(result.points.size()>=config.maximum_points)break;
    }
    return result;
}

ComplexMatrix single_ended_to_mixed_mode(const ComplexMatrix& s) {
    if (s.size() != 4U) throw std::invalid_argument("mixed-mode conversion requires a 4-port network");
    const double inv_root2 = 1.0 / std::sqrt(2.0);
    ComplexMatrix c(4U);
    c(0,0)=inv_root2; c(0,1)=-inv_root2;
    c(1,2)=inv_root2; c(1,3)=-inv_root2;
    c(2,0)=inv_root2; c(2,1)= inv_root2;
    c(3,2)=inv_root2; c(3,3)= inv_root2;
    // C is real orthogonal and symmetric only in pair sub-blocks; use transpose explicitly.
    ComplexMatrix ct(4U);
    for (std::size_t r=0;r<4U;++r) for(std::size_t col=0;col<4U;++col) ct(r,col)=c(col,r);
    return multiply(multiply(c,s),ct);
}

LMatch synthesize_l_match(double rs, double rl) {
    if (!(rs > 0.0) || !(rl > 0.0) || !std::isfinite(rs) || !std::isfinite(rl))
        throw std::invalid_argument("L-match resistances must be finite and positive");
    if (std::abs(rs - rl) < 1.0e-15 * std::max(rs, rl)) return {};
    if (rl > rs) {
        const double q = std::sqrt(rl / rs - 1.0);
        return {rs * q, q / rl, false};
    }
    const double q = std::sqrt(rs / rl - 1.0);
    // High-to-low low-pass dual uses the shunt element at the source side.
    return {rl * q, q / rs, true};
}

} // namespace cfd::rf

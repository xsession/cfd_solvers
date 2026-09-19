#include "cfd/rf/modal_elements.hpp"

#include <cmath>
#include <stdexcept>

namespace cfd::rf {

Matrix2C series_modal_discontinuity(Complex impedance) {
    if (!std::isfinite(impedance.real()) || !std::isfinite(impedance.imag()))
        throw std::invalid_argument("invalid modal series impedance");
    return {1.0, impedance, 0.0, 1.0};
}

Matrix2C shunt_modal_discontinuity(Complex admittance) {
    if (!std::isfinite(admittance.real()) || !std::isfinite(admittance.imag()))
        throw std::invalid_argument("invalid modal shunt admittance");
    return {1.0, 0.0, admittance, 1.0};
}

Matrix2C modal_discontinuity_s(Complex series_impedance, Complex shunt_admittance, double reference_impedance) {
    if (!(reference_impedance > 0.0) || !std::isfinite(reference_impedance))
        throw std::invalid_argument("invalid modal reference impedance");
    return abcd_to_s(
        cascade_abcd(series_modal_discontinuity(series_impedance), shunt_modal_discontinuity(shunt_admittance)),
        reference_impedance);
}

} // namespace cfd::rf

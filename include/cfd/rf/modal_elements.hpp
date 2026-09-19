#pragma once

#include "cfd/rf/network.hpp"

namespace cfd::rf {

// Equivalent two-port modal discontinuities. These ABCD elements can be
// cascaded with field-extracted ports or converted to S parameters.
[[nodiscard]] Matrix2C series_modal_discontinuity(Complex impedance_ohm);
[[nodiscard]] Matrix2C shunt_modal_discontinuity(Complex admittance_siemens);
[[nodiscard]] Matrix2C modal_discontinuity_s(Complex series_impedance_ohm, Complex shunt_admittance_siemens,
                                             double reference_impedance = 50.0);

} // namespace cfd::rf

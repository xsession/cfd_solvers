#pragma once

#include "cfd/rf/nport.hpp"

#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

namespace cfd::rf {

struct ThinWireMomConfig {
    double length_m{};
    double radius_m{};
    double frequency_hz{};
    std::size_t segments{31U}; // odd: feed is the center segment
    double feed_voltage_v{1.0};
    std::size_t quadrature_order{8U};
};

struct ThinWireMomResult {
    std::vector<double> segment_center_z_m;
    std::vector<std::complex<double>> current_a;
    std::complex<double> feed_impedance_ohm{};
    double feed_current_a{};
};

// Thin straight PEC wire along z, pulse basis + point matching applied to
// Pocklington's thin-wire EFIE. A delta-gap electric field is imposed on the
// center segment. This is a readable reference solver, not a NEC replacement.
[[nodiscard]] ThinWireMomResult solve_center_fed_thin_wire(const ThinWireMomConfig& config);

struct ParallelThinWire {
    double x_m{};
    double y_m{};
    double center_z_m{};
    double length_m{};
    double radius_m{};
    std::size_t segments{31U};
};

struct ThinWireDeltaGapFeed {
    std::size_t wire{};
    // max() selects the center segment and therefore requires odd segmentation.
    std::size_t segment{std::numeric_limits<std::size_t>::max()};
    std::complex<double> voltage_v{1.0,0.0};
};

struct ThinWireSeriesLoad {
    std::size_t wire{};
    std::size_t segment{};
    std::complex<double> impedance_ohm{};
};

// A two-port transmission-line section coupled between two MoM segment-current
// unknowns. Port-current signs map the local segment-current orientation into
// current entering each transmission-line port and must be +1 or -1.
struct ThinWireTransmissionLineLoad {
    std::size_t first_wire{};
    std::size_t first_segment{};
    std::size_t second_wire{};
    std::size_t second_segment{};
    std::complex<double> characteristic_impedance_ohm{50.0,0.0};
    std::complex<double> propagation_constant_per_m{}; // alpha + j*beta
    double length_m{};
    int first_current_sign{1};
    int second_current_sign{1};
};

struct ParallelWireMomConfig {
    double frequency_hz{};
    std::vector<ParallelThinWire> wires;
    std::vector<ThinWireDeltaGapFeed> feeds;
    std::vector<ThinWireSeriesLoad> loads;
    std::size_t quadrature_order{8U};
};

struct WirePoint3 {
    double x{};
    double y{};
    double z{};
};

struct OrientedThinWire {
    WirePoint3 start_m;
    WirePoint3 end_m;
    double radius_m{};
    std::size_t segments{31U};
    // Zero conductivity means PEC. Positive conductivity adds the high-
    // frequency round-wire surface-resistance baseline.
    double conductivity_s_per_m{};
    double relative_permeability{1.0};
    // Optional equivalent distributed dielectric loss in ohm/m. This is a
    // circuit-equivalent loss term, not a replacement for volumetric dielectric fields.
    double dielectric_loss_ohm_per_m{};
};

struct OrientedWireSegment {
    std::size_t wire{};
    std::size_t local_segment{};
    WirePoint3 center_m;
    WirePoint3 tangent;
    double length_m{};
};

enum class ThinWireGroundModel {
    free_space,
    pec_plane
};

struct OrientedWireMomConfig {
    double frequency_hz{};
    std::vector<OrientedThinWire> wires;
    std::vector<ThinWireDeltaGapFeed> feeds;
    std::vector<ThinWireSeriesLoad> loads;
    std::vector<ThinWireTransmissionLineLoad> transmission_line_loads;
    std::size_t quadrature_order{8U};
    // Zero preserves the legacy disjoint-wire rule. A positive tolerance enables
    // endpoint-only junctions; crossing/interior intersections remain rejected.
    double junction_tolerance_m{};
    ThinWireGroundModel ground_model{ThinWireGroundModel::free_space};
    double ground_plane_z_m{};
};

struct OrientedWireMomResult {
    std::vector<OrientedWireSegment> segments;
    std::vector<std::complex<double>> current_a;
    std::vector<std::complex<double>> feed_impedance_ohm;
    std::vector<double> feed_current_a;
    double accepted_power_w{};
    double conductor_loss_w{};
    double dielectric_loss_w{};
    double lumped_load_loss_w{};
    double transmission_line_loss_w{};
    double radiated_power_w{};
    double radiation_efficiency{};
};

// Arbitrarily oriented, mutually coupled straight PEC wires. The dyadic free-
// space Green kernel is projected onto source and observation tangents. Wires
// may be joined exactly at endpoints when junction_tolerance_m is positive.
// Junction KCL is imposed with a constrained saddle-point solve; interior crossings remain rejected.
[[nodiscard]] OrientedWireMomResult solve_oriented_thin_wires(const OrientedWireMomConfig& config);

// Multiport impedance extraction using the configured delta-gap feed locations.
// Each port is excited independently with 1 V while all other impressed feed
// voltages are zero; the resulting port-current admittance matrix is inverted.
[[nodiscard]] ComplexMatrix oriented_wire_port_impedance_matrix(const OrientedWireMomConfig& config);
[[nodiscard]] ComplexMatrix oriented_wire_port_s_parameters(const OrientedWireMomConfig& config,
                                                            std::span<const double> reference_impedance);

struct ParallelWireSegment {
    std::size_t wire{};
    std::size_t local_segment{};
    double x_m{};
    double y_m{};
    double center_z_m{};
    double length_m{};
};

struct ParallelWireMomResult {
    std::vector<ParallelWireSegment> segments;
    std::vector<std::complex<double>> current_a;
    std::vector<std::complex<double>> feed_impedance_ohm;
    std::vector<double> feed_current_a;
};

// Compatibility wrapper for the z-directed parallel-wire API. It delegates to
// the arbitrary-orientation solver while preserving the original result layout.
[[nodiscard]] ParallelWireMomResult solve_parallel_thin_wires(const ParallelWireMomConfig& config);

} // namespace cfd::rf

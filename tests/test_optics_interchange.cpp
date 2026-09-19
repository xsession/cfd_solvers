#include "cfd/optics/interchange.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

cfd::optics::SequentialOpticalSystem sample_system() {
    using namespace cfd::optics;
    SequentialOpticalSystem system;
    system.add_surface({SurfaceType::plane, 0.0, 0.0, 0.25, 1.5, 0.0, {}});
    system.add_surface({SurfaceType::sphere, 0.02, 0.10, 0.20, 1.0, 0.0, {}});
    system.add_surface({SurfaceType::even_asphere, 0.05, -0.12, 0.18, 1.0, -1.0, {1.0e-5, -2.0e-7, 0.0, 0.0}});
    return system;
}

void round_trip(cfd::optics::LensInterchangeFormat format) {
    using namespace cfd::optics;
    const auto source = sample_system();
    const auto text = export_lens_system(source, format);
    const auto restored = import_lens_system(text, format);
    require(restored.surfaces().size() == source.surfaces().size(), "lens surface count round trip");
    for (std::size_t i = 0; i < source.surfaces().size(); ++i) {
        const auto& a = source.surfaces()[i];
        const auto& b = restored.surfaces()[i];
        require(a.type == b.type, "lens surface type round trip");
        require(std::abs(a.vertex_z - b.vertex_z) < 1.0e-12, "lens vertex position round trip");
        require(std::abs(a.aperture_radius - b.aperture_radius) < 1.0e-12, "lens aperture round trip");
        require(std::abs(a.refractive_index_after - b.refractive_index_after) < 1.0e-12, "lens material round trip");
    }
}

} // namespace

int main() {
    try {
        round_trip(cfd::optics::LensInterchangeFormat::zemax_zmx);
        round_trip(cfd::optics::LensInterchangeFormat::codev_seq);
        round_trip(cfd::optics::LensInterchangeFormat::oslo_len);
        std::cout << "optics lens interchange regression passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "optics interchange regression failed: " << error.what() << '\n';
        return 1;
    }
}

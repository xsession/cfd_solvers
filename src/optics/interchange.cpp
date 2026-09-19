#include "cfd/optics/interchange.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfd::optics {
namespace {

struct SurfaceRecord {
    bool active{};
    SequentialSurface surface{};
    double thickness{};
    double curvature{};
};

[[nodiscard]] const char* type_word(SurfaceType type, LensInterchangeFormat format) {
    if (type == SurfaceType::plane)
        return format == LensInterchangeFormat::codev_seq ? "PLANO" : "PLANE";
    if (type == SurfaceType::sphere)
        return format == LensInterchangeFormat::codev_seq ? "STANDARD" : "STANDARD";
    if (type == SurfaceType::conic)
        return format == LensInterchangeFormat::codev_seq ? "CONIC" : "CONIC";
    return format == LensInterchangeFormat::codev_seq ? "EVENASPH" : "EVENASPH";
}

[[nodiscard]] double parse_last_number(const std::vector<std::string>& tokens) {
    for (auto it = tokens.rbegin(); it != tokens.rend(); ++it) {
        try {
            std::size_t used = 0U;
            const double value = std::stod(*it, &used);
            if (used == it->size())
                return value;
        } catch (const std::exception&) {
        }
    }
    throw std::invalid_argument("lens record is missing a numeric value");
}

void finish_surface(std::vector<SurfaceRecord>& records, SurfaceRecord& record, double& previous_z) {
    if (!record.active)
        return;
    if (record.surface.type != SurfaceType::plane && std::abs(record.curvature) > 0.0)
        record.surface.radius = 1.0 / record.curvature;
    if (record.surface.type != SurfaceType::plane && record.surface.radius == 0.0)
        record.surface.type = SurfaceType::plane;
    record.surface.vertex_z = previous_z + record.thickness;
    previous_z = record.surface.vertex_z;
    records.push_back(record);
    record = {};
}

} // namespace

std::string export_lens_system(const SequentialOpticalSystem& system, LensInterchangeFormat format) {
    std::ostringstream output;
    output << std::setprecision(17) << "# cfd_solvers sequential lens interchange\n";
    output << "# format "
           << (format == LensInterchangeFormat::zemax_zmx   ? "ZMX"
               : format == LensInterchangeFormat::codev_seq ? "SEQ"
                                                            : "LEN")
           << "\n";
    double previous_z = 0.0;
    for (std::size_t index = 0; index < system.surfaces().size(); ++index) {
        const auto& surface = system.surfaces()[index];
        const double thickness = surface.vertex_z - previous_z;
        const double curvature = surface.type == SurfaceType::plane ? 0.0 : 1.0 / surface.radius;
        output << "SURF " << index << "\n";
        output << "TYPE " << type_word(surface.type, format) << "\n";
        if (format == LensInterchangeFormat::zemax_zmx) {
            output << "CURV " << curvature << "\nDISZ " << thickness << "\nDIAM " << surface.aperture_radius
                   << "\nGLAS " << surface.refractive_index_after << "\nCONI " << surface.conic_constant << "\n";
            for (std::size_t i = 0; i < surface.even_coefficients.size(); ++i)
                output << "PARM " << (i + 1U) << ' ' << surface.even_coefficients[i] << "\n";
        } else if (format == LensInterchangeFormat::codev_seq) {
            output << "RADIUS " << (surface.type == SurfaceType::plane ? 0.0 : surface.radius) << "\nTHICK "
                   << thickness << "\nAPERTURE " << surface.aperture_radius << "\nGLASS "
                   << surface.refractive_index_after << "\nCONIC " << surface.conic_constant << "\n";
            for (std::size_t i = 0; i < surface.even_coefficients.size(); ++i)
                output << "ASPH " << (i + 1U) << ' ' << surface.even_coefficients[i] << "\n";
        } else {
            output << "RD " << (surface.type == SurfaceType::plane ? 0.0 : surface.radius) << "\nTH " << thickness
                   << "\nAP " << surface.aperture_radius << "\nGLASS " << surface.refractive_index_after << "\nCONI "
                   << surface.conic_constant << "\n";
            for (std::size_t i = 0; i < surface.even_coefficients.size(); ++i)
                output << "COEF " << (i + 1U) << ' ' << surface.even_coefficients[i] << "\n";
        }
        previous_z = surface.vertex_z;
    }
    return output.str();
}

SequentialOpticalSystem import_lens_system(std::string_view text, LensInterchangeFormat format) {
    std::istringstream input{std::string(text)};
    std::string line;
    std::vector<SurfaceRecord> records;
    SurfaceRecord record;
    double previous_z = 0.0;
    while (std::getline(input, line)) {
        const auto comment = line.find_first_of("#!");
        if (comment != std::string::npos)
            line.resize(comment);
        std::istringstream tokens_stream(line);
        std::vector<std::string> tokens;
        std::string token;
        while (tokens_stream >> token)
            tokens.push_back(token);
        if (tokens.empty())
            continue;
        std::string keyword = tokens.front();
        std::transform(keyword.begin(), keyword.end(), keyword.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (keyword == "SURF" || keyword == "SURFACE") {
            finish_surface(records, record, previous_z);
            record.active = true;
            continue;
        }
        if (!record.active)
            continue;
        if (keyword == "TYPE" && tokens.size() >= 2U) {
            std::string type = tokens[1];
            std::transform(type.begin(), type.end(), type.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (type == "PLANE" || type == "PLANO")
                record.surface.type = SurfaceType::plane;
            else if (type == "CONIC")
                record.surface.type = SurfaceType::conic;
            else if (type == "EVENASPH" || type == "EVEN_ASPHERE")
                record.surface.type = SurfaceType::even_asphere;
            else
                record.surface.type = SurfaceType::sphere;
        } else if (keyword == "CURV") {
            record.curvature = parse_last_number(tokens);
        } else if (keyword == "RADIUS" || keyword == "RD") {
            const double radius = parse_last_number(tokens);
            record.curvature = radius == 0.0 ? 0.0 : 1.0 / radius;
        } else if (keyword == "DISZ" || keyword == "THICK" || keyword == "TH") {
            record.thickness = parse_last_number(tokens);
        } else if (keyword == "DIAM" || keyword == "APERTURE" || keyword == "AP") {
            record.surface.aperture_radius = parse_last_number(tokens);
        } else if (keyword == "GLAS" || keyword == "GLASS") {
            record.surface.refractive_index_after = parse_last_number(tokens);
        } else if (keyword == "CONI") {
            record.surface.conic_constant = parse_last_number(tokens);
        } else if (keyword == "PARM" || keyword == "ASPH" || keyword == "COEF") {
            if (tokens.size() >= 3U) {
                const double index = parse_last_number({tokens[1]});
                const double value = parse_last_number({tokens[2]});
                if (index >= 1.0 && index <= 4.0)
                    record.surface.even_coefficients[static_cast<std::size_t>(index - 1.0)] = value;
            }
        }
    }
    finish_surface(records, record, previous_z);
    if (records.empty())
        throw std::invalid_argument("lens interchange contains no surfaces");
    SequentialOpticalSystem system;
    for (const auto& item : records)
        system.add_surface(item.surface);
    (void)format;
    return system;
}

} // namespace cfd::optics

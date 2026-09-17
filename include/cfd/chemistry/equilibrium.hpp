#pragma once
namespace cfd::chemistry {
struct BinaryPrecipitationResult {
    double dissolved_a{};
    double dissolved_b{};
    double solid_amount{};
    double reaction_extent{};
    double saturation_ratio{};
};
[[nodiscard]] BinaryPrecipitationResult equilibrate_binary_salt(double dissolved_a,double dissolved_b,double solid_amount,double solubility_product);
}

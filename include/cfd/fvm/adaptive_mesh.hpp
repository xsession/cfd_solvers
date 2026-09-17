#pragma once

#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace cfd::fvm {

struct AdaptiveHexCell {
    Vec3 lower{};
    Vec3 upper{};
    std::size_t level{};
    std::size_t root{};
    std::uint64_t lineage{};
};

class AdaptiveHexMesh {
public:
    AdaptiveHexMesh(std::vector<AdaptiveHexCell> cells, Vec3 domain_lower, Vec3 domain_upper);

    [[nodiscard]] static AdaptiveHexMesh cartesian(std::size_t nx,
                                                   std::size_t ny,
                                                   std::size_t nz,
                                                   double lx = 1.0,
                                                   double ly = 1.0,
                                                   double lz = 1.0);

    [[nodiscard]] const std::vector<AdaptiveHexCell>& cells() const noexcept { return cells_; }
    [[nodiscard]] std::size_t cell_count() const noexcept { return cells_.size(); }
    [[nodiscard]] Vec3 domain_lower() const noexcept { return domain_lower_; }
    [[nodiscard]] Vec3 domain_upper() const noexcept { return domain_upper_; }
    [[nodiscard]] PolyMesh poly_mesh() const;
    [[nodiscard]] double total_volume() const noexcept;

private:
    std::vector<AdaptiveHexCell> cells_;
    Vec3 domain_lower_{};
    Vec3 domain_upper_{};

    void validate() const;
};

struct CellOverlapWeight {
    std::size_t source_cell{};
    double volume{};
};

struct TopologyChangeMap {
    std::size_t old_cell_count{};
    std::size_t new_cell_count{};
    std::vector<std::vector<CellOverlapWeight>> new_from_old;
};

struct AdaptiveTopologyResult {
    AdaptiveHexMesh mesh;
    TopologyChangeMap mapping;
};

class MeshTopologyOperation {
public:
    virtual ~MeshTopologyOperation() = default;
    [[nodiscard]] virtual AdaptiveTopologyResult apply(const AdaptiveHexMesh& mesh) const = 0;
};

class MarkedHexRefinement final : public MeshTopologyOperation {
public:
    MarkedHexRefinement(std::vector<unsigned char> marked,
                        std::size_t maximum_level,
                        bool enforce_two_to_one = true);
    [[nodiscard]] AdaptiveTopologyResult apply(const AdaptiveHexMesh& mesh) const override;

private:
    std::vector<unsigned char> marked_;
    std::size_t maximum_level_{};
    bool enforce_two_to_one_{true};
};

class MarkedHexCoarsening final : public MeshTopologyOperation {
public:
    MarkedHexCoarsening(std::vector<unsigned char> marked,
                        std::size_t minimum_level = 0U,
                        bool enforce_two_to_one = true);
    [[nodiscard]] AdaptiveTopologyResult apply(const AdaptiveHexMesh& mesh) const override;

private:
    std::vector<unsigned char> marked_;
    std::size_t minimum_level_{};
    bool enforce_two_to_one_{true};
};

[[nodiscard]] TopologyChangeMap build_topology_change_map(const AdaptiveHexMesh& old_mesh,
                                                           const AdaptiveHexMesh& new_mesh);

[[nodiscard]] std::vector<double> conservative_adaptive_remap(const AdaptiveHexMesh& old_mesh,
                                                              std::span<const double> old_field,
                                                              const AdaptiveHexMesh& new_mesh);

[[nodiscard]] std::vector<double> scalar_jump_error_indicator(const AdaptiveHexMesh& mesh,
                                                              std::span<const double> field);

[[nodiscard]] std::vector<unsigned char> mark_dorfler(std::span<const double> indicators,
                                                      double theta);

struct AdaptiveAmrConfig {
    std::size_t minimum_level{0U};
    std::size_t maximum_level{4U};
    double refine_theta{0.5};
    double coarsen_relative_threshold{0.05};
    bool enforce_two_to_one{true};
};

struct AdaptiveAmrResult {
    AdaptiveHexMesh mesh;
    std::vector<double> field;
    std::vector<double> indicator;
    std::size_t refined_cell_count{};
    std::size_t coarsened_cell_count{};
};

[[nodiscard]] AdaptiveAmrResult adapt_scalar_field(const AdaptiveHexMesh& mesh,
                                                   std::span<const double> field,
                                                   const AdaptiveAmrConfig& config = {});

[[nodiscard]] bool satisfies_two_to_one(const AdaptiveHexMesh& mesh);

} // namespace cfd::fvm

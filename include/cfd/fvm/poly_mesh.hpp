#pragma once

#include <cstddef>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::fvm {

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

[[nodiscard]] constexpr Vec3 operator+(Vec3 a, Vec3 b) noexcept { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
[[nodiscard]] constexpr Vec3 operator-(Vec3 a, Vec3 b) noexcept { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
[[nodiscard]] constexpr Vec3 operator*(Vec3 a, double s) noexcept { return {a.x*s,a.y*s,a.z*s}; }
[[nodiscard]] constexpr Vec3 operator*(double s, Vec3 a) noexcept { return a*s; }
[[nodiscard]] constexpr Vec3 operator/(Vec3 a, double s) noexcept { return {a.x/s,a.y/s,a.z/s}; }
constexpr Vec3& operator+=(Vec3& a, Vec3 b) noexcept { a=a+b; return a; }
constexpr Vec3& operator-=(Vec3& a, Vec3 b) noexcept { a=a-b; return a; }
[[nodiscard]] constexpr double dot(Vec3 a, Vec3 b) noexcept { return a.x*b.x+a.y*b.y+a.z*b.z; }
[[nodiscard]] double magnitude(Vec3 a) noexcept;

inline constexpr std::size_t invalid_cell = std::numeric_limits<std::size_t>::max();
inline constexpr std::size_t invalid_patch = std::numeric_limits<std::size_t>::max();

struct Cell {
    Vec3 center{};
    double volume{};
};

struct Face {
    std::size_t owner{};
    std::size_t neighbour{invalid_cell};
    Vec3 center{};
    Vec3 area{}; // points out of owner; for internal faces it points owner -> neighbour
    std::size_t patch{invalid_patch};

    [[nodiscard]] bool boundary() const noexcept { return neighbour == invalid_cell; }
};

struct BoundaryPatch {
    std::string name;
};

class PolyMesh {
public:
    PolyMesh(std::vector<Cell> cells,
             std::vector<Face> faces,
             std::vector<BoundaryPatch> patches);

    [[nodiscard]] const std::vector<Cell>& cells() const noexcept { return cells_; }
    [[nodiscard]] const std::vector<Face>& faces() const noexcept { return faces_; }
    [[nodiscard]] const std::vector<BoundaryPatch>& patches() const noexcept { return patches_; }
    [[nodiscard]] const std::vector<std::vector<std::size_t>>& cell_faces() const noexcept { return cell_faces_; }
    [[nodiscard]] std::size_t cell_count() const noexcept { return cells_.size(); }
    [[nodiscard]] std::size_t face_count() const noexcept { return faces_.size(); }
    [[nodiscard]] std::size_t boundary_face_count() const noexcept;
    [[nodiscard]] std::size_t patch_index(std::string_view name) const;

private:
    std::vector<Cell> cells_;
    std::vector<Face> faces_;
    std::vector<BoundaryPatch> patches_;
    std::vector<std::vector<std::size_t>> cell_faces_;

    void validate_and_build_adjacency();
};

[[nodiscard]] PolyMesh make_cartesian_hexa_mesh(std::size_t nx,
                                                std::size_t ny,
                                                std::size_t nz,
                                                double lx = 1.0,
                                                double ly = 1.0,
                                                double lz = 1.0);

// Affine xy shear of the Cartesian mesh. Useful for deterministic
// non-orthogonal finite-volume regression tests.
[[nodiscard]] PolyMesh make_sheared_cartesian_hexa_mesh(std::size_t nx,
                                                        std::size_t ny,
                                                        std::size_t nz,
                                                        double lx = 1.0,
                                                        double ly = 1.0,
                                                        double lz = 1.0,
                                                        double shear_xy = 0.25);

} // namespace cfd::fvm

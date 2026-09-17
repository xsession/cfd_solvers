#include "cfd/io/fvm_checkpoint.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <regex>
#include <stdexcept>
#include <utility>

#ifdef CFD_HAS_HDF5
#include <hdf5.h>
#endif

namespace cfd::io {
namespace {
void validate_field_name(const std::string& name) {
    static const std::regex valid("[A-Za-z0-9_.-]+");
    if (name.empty() || !std::regex_match(name, valid))
        throw std::invalid_argument("FVM HDF5 field names are restricted to [A-Za-z0-9_.-]+");
}

void validate_inputs(const cfd::fvm::PolyMesh& mesh,
                     std::span<const FvmScalarFieldView> fields,
                     double time) {
    if (!std::isfinite(time) || time < 0.0) throw std::invalid_argument("invalid FVM checkpoint time");
    std::vector<std::string> names;
    names.reserve(fields.size());
    for (const auto& field : fields) {
        validate_field_name(field.name);
        if (field.values.size() != mesh.cell_count())
            throw std::invalid_argument("FVM checkpoint scalar field size mismatch");
        if (!std::all_of(field.values.begin(), field.values.end(), [](double value) { return std::isfinite(value); }))
            throw std::invalid_argument("FVM checkpoint scalar field contains a non-finite value");
        if (std::find(names.begin(), names.end(), field.name) != names.end())
            throw std::invalid_argument("duplicate FVM checkpoint field name");
        names.push_back(field.name);
    }
}

#ifndef CFD_HAS_HDF5
[[noreturn]] void unavailable() {
    throw std::runtime_error("FVM HDF5 checkpoint I/O requires CFD_ENABLE_HDF5 and an HDF5 library");
}
#else
class H5Handle {
public:
    H5Handle() = default;
    H5Handle(hid_t value, herr_t (*closer)(hid_t)) : value_(value), closer_(closer) {}
    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;
    H5Handle(H5Handle&& other) noexcept : value_(std::exchange(other.value_, -1)), closer_(other.closer_) {}
    H5Handle& operator=(H5Handle&& other) noexcept {
        if (this != &other) {
            reset();
            value_ = std::exchange(other.value_, -1);
            closer_ = other.closer_;
        }
        return *this;
    }
    ~H5Handle() { reset(); }
    [[nodiscard]] hid_t get() const noexcept { return value_; }
    explicit operator bool() const noexcept { return value_ >= 0; }
private:
    void reset() noexcept { if (value_ >= 0 && closer_) closer_(value_); value_ = -1; }
    hid_t value_{-1};
    herr_t (*closer_)(hid_t){};
};

void require_h5(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

void create_group(hid_t file, const char* path) {
    H5Handle group(H5Gcreate2(file, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    require_h5(static_cast<bool>(group), "cannot create FVM HDF5 group");
}

void write_double_dataset(hid_t file, const std::string& path, const double* data, std::span<const hsize_t> dims) {
    H5Handle space(H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr), H5Sclose);
    require_h5(static_cast<bool>(space), "cannot create FVM HDF5 dataspace");
    H5Handle set(H5Dcreate2(file, path.c_str(), H5T_IEEE_F64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "cannot create FVM HDF5 dataset");
    require_h5(H5Dwrite(set.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) >= 0,
               "cannot write FVM HDF5 double dataset");
}

void write_u64_dataset(hid_t file, const std::string& path, const std::uint64_t* data, std::span<const hsize_t> dims) {
    H5Handle space(H5Screate_simple(static_cast<int>(dims.size()), dims.data(), nullptr), H5Sclose);
    require_h5(static_cast<bool>(space), "cannot create FVM HDF5 dataspace");
    H5Handle set(H5Dcreate2(file, path.c_str(), H5T_STD_U64LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "cannot create FVM HDF5 dataset");
    require_h5(H5Dwrite(set.get(), H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) >= 0,
               "cannot write FVM HDF5 integer dataset");
}

void write_string_table(hid_t file, const std::string& path, std::span<const std::string> strings) {
    const std::size_t width = std::max<std::size_t>(1U, [&] {
        std::size_t out = 1U;
        for (const auto& value : strings) out = std::max(out, value.size() + 1U);
        return out;
    }());
    std::vector<char> bytes(strings.size() * width, '\0');
    for (std::size_t i = 0; i < strings.size(); ++i)
        std::memcpy(bytes.data() + i * width, strings[i].data(), strings[i].size());
    const hsize_t dims[2]{static_cast<hsize_t>(strings.size()), static_cast<hsize_t>(width)};
    H5Handle space(H5Screate_simple(2, dims, nullptr), H5Sclose);
    require_h5(static_cast<bool>(space), "cannot create FVM HDF5 string dataspace");
    H5Handle set(H5Dcreate2(file, path.c_str(), H5T_STD_I8LE, space.get(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "cannot create FVM HDF5 string dataset");
    require_h5(H5Dwrite(set.get(), H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data()) >= 0,
               "cannot write FVM HDF5 string dataset");
}

std::vector<hsize_t> dataset_dims(hid_t dataset) {
    H5Handle space(H5Dget_space(dataset), H5Sclose);
    require_h5(static_cast<bool>(space), "cannot inspect FVM HDF5 dataset");
    const int rank = H5Sget_simple_extent_ndims(space.get());
    require_h5(rank >= 0, "cannot inspect FVM HDF5 dataset rank");
    std::vector<hsize_t> dims(static_cast<std::size_t>(rank));
    require_h5(H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr) >= 0,
               "cannot inspect FVM HDF5 dataset dimensions");
    return dims;
}

std::vector<double> read_double_dataset(hid_t file, const std::string& path, std::span<const hsize_t> expected) {
    H5Handle set(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "missing FVM HDF5 double dataset");
    const auto dims = dataset_dims(set.get());
    require_h5(dims.size() == expected.size() && std::equal(dims.begin(), dims.end(), expected.begin()),
               "unexpected FVM HDF5 double dataset shape");
    std::size_t count = 1U;
    for (const auto dim : dims) count *= static_cast<std::size_t>(dim);
    std::vector<double> values(count);
    require_h5(H5Dread(set.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0,
               "cannot read FVM HDF5 double dataset");
    return values;
}

std::vector<std::uint64_t> read_u64_dataset(hid_t file, const std::string& path, std::span<const hsize_t> expected) {
    H5Handle set(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "missing FVM HDF5 integer dataset");
    const auto dims = dataset_dims(set.get());
    require_h5(dims.size() == expected.size() && std::equal(dims.begin(), dims.end(), expected.begin()),
               "unexpected FVM HDF5 integer dataset shape");
    std::size_t count = 1U;
    for (const auto dim : dims) count *= static_cast<std::size_t>(dim);
    std::vector<std::uint64_t> values(count);
    require_h5(H5Dread(set.get(), H5T_NATIVE_UINT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0,
               "cannot read FVM HDF5 integer dataset");
    return values;
}

std::vector<std::string> read_string_table(hid_t file, const std::string& path) {
    H5Handle set(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(set), "missing FVM HDF5 string table");
    const auto dims = dataset_dims(set.get());
    require_h5(dims.size() == 2U && dims[1] > 0U, "unexpected FVM HDF5 string-table shape");
    const std::size_t rows = static_cast<std::size_t>(dims[0]);
    const std::size_t width = static_cast<std::size_t>(dims[1]);
    std::vector<char> bytes(rows * width, '\0');
    require_h5(H5Dread(set.get(), H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data()) >= 0,
               "cannot read FVM HDF5 string table");
    std::vector<std::string> strings;
    strings.reserve(rows);
    for (std::size_t row = 0; row < rows; ++row) {
        const char* begin = bytes.data() + row * width;
        const char* end = std::find(begin, begin + width, '\0');
        strings.emplace_back(begin, static_cast<std::size_t>(end - begin));
    }
    return strings;
}
#endif
} // namespace

bool fvm_hdf5_checkpoint_available() noexcept {
#ifdef CFD_HAS_HDF5
    return true;
#else
    return false;
#endif
}

void write_fvm_hdf5_checkpoint(const std::filesystem::path& path,
                               const cfd::fvm::PolyMesh& mesh,
                               std::span<const FvmScalarFieldView> scalar_fields,
                               double time,
                               std::size_t step) {
    validate_inputs(mesh, scalar_fields, time);
#ifndef CFD_HAS_HDF5
    static_cast<void>(path); static_cast<void>(step); unavailable();
#else
    H5Handle file(H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    require_h5(static_cast<bool>(file), "cannot create FVM HDF5 checkpoint");
    create_group(file.get(), "/Mesh");
    create_group(file.get(), "/Fields");
    create_group(file.get(), "/Meta");

    std::vector<double> cells;
    cells.reserve(mesh.cell_count() * 4U);
    for (const auto& cell : mesh.cells()) {
        cells.insert(cells.end(), {cell.center.x, cell.center.y, cell.center.z, cell.volume});
    }
    const hsize_t cell_dims[2]{static_cast<hsize_t>(mesh.cell_count()), 4U};
    write_double_dataset(file.get(), "/Mesh/Cells", cells.data(), cell_dims);

    std::vector<double> face_geometry;
    std::vector<std::uint64_t> face_index;
    face_geometry.reserve(mesh.face_count() * 6U);
    face_index.reserve(mesh.face_count() * 3U);
    constexpr std::uint64_t invalid = std::numeric_limits<std::uint64_t>::max();
    for (const auto& face : mesh.faces()) {
        face_geometry.insert(face_geometry.end(), {face.center.x, face.center.y, face.center.z,
                                                   face.area.x, face.area.y, face.area.z});
        face_index.push_back(static_cast<std::uint64_t>(face.owner));
        face_index.push_back(face.boundary() ? invalid : static_cast<std::uint64_t>(face.neighbour));
        face_index.push_back(face.patch == cfd::fvm::invalid_patch ? invalid : static_cast<std::uint64_t>(face.patch));
    }
    const hsize_t face_geom_dims[2]{static_cast<hsize_t>(mesh.face_count()), 6U};
    const hsize_t face_index_dims[2]{static_cast<hsize_t>(mesh.face_count()), 3U};
    write_double_dataset(file.get(), "/Mesh/FaceGeometry", face_geometry.data(), face_geom_dims);
    write_u64_dataset(file.get(), "/Mesh/FaceIndex", face_index.data(), face_index_dims);

    std::vector<std::string> patch_names;
    patch_names.reserve(mesh.patches().size());
    for (const auto& patch : mesh.patches()) patch_names.push_back(patch.name);
    write_string_table(file.get(), "/Mesh/PatchNames", patch_names);

    const hsize_t one[1]{1U};
    write_double_dataset(file.get(), "/Meta/Time", &time, one);
    const std::uint64_t step64 = static_cast<std::uint64_t>(step);
    write_u64_dataset(file.get(), "/Meta/Step", &step64, one);
    std::vector<std::string> field_names;
    field_names.reserve(scalar_fields.size());
    for (const auto& field : scalar_fields) field_names.push_back(field.name);
    write_string_table(file.get(), "/Meta/FieldNames", field_names);
    const hsize_t field_dims[1]{static_cast<hsize_t>(mesh.cell_count())};
    for (const auto& field : scalar_fields)
        write_double_dataset(file.get(), "/Fields/" + field.name, field.values.data(), field_dims);
#endif
}

FvmCheckpoint read_fvm_hdf5_checkpoint(const std::filesystem::path& path) {
#ifndef CFD_HAS_HDF5
    static_cast<void>(path); unavailable();
#else
    H5Handle file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    require_h5(static_cast<bool>(file), "cannot open FVM HDF5 checkpoint");

    H5Handle cells_set(H5Dopen2(file.get(), "/Mesh/Cells", H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(cells_set), "missing FVM checkpoint cells");
    const auto cell_shape = dataset_dims(cells_set.get());
    require_h5(cell_shape.size() == 2U && cell_shape[1] == 4U, "invalid FVM checkpoint cell shape");
    const std::size_t cell_count = static_cast<std::size_t>(cell_shape[0]);
    const hsize_t cell_dims[2]{cell_shape[0], 4U};
    const auto cell_data = read_double_dataset(file.get(), "/Mesh/Cells", cell_dims);
    std::vector<cfd::fvm::Cell> cells;
    cells.reserve(cell_count);
    for (std::size_t i = 0; i < cell_count; ++i)
        cells.push_back({{cell_data[4U*i], cell_data[4U*i+1U], cell_data[4U*i+2U]}, cell_data[4U*i+3U]});

    H5Handle faces_set(H5Dopen2(file.get(), "/Mesh/FaceGeometry", H5P_DEFAULT), H5Dclose);
    require_h5(static_cast<bool>(faces_set), "missing FVM checkpoint faces");
    const auto face_shape = dataset_dims(faces_set.get());
    require_h5(face_shape.size() == 2U && face_shape[1] == 6U, "invalid FVM checkpoint face shape");
    const std::size_t face_count = static_cast<std::size_t>(face_shape[0]);
    const hsize_t face_geom_dims[2]{face_shape[0], 6U};
    const hsize_t face_index_dims[2]{face_shape[0], 3U};
    const auto face_geometry = read_double_dataset(file.get(), "/Mesh/FaceGeometry", face_geom_dims);
    const auto face_index = read_u64_dataset(file.get(), "/Mesh/FaceIndex", face_index_dims);
    constexpr std::uint64_t invalid = std::numeric_limits<std::uint64_t>::max();
    std::vector<cfd::fvm::Face> faces;
    faces.reserve(face_count);
    for (std::size_t i = 0; i < face_count; ++i) {
        const auto owner = static_cast<std::size_t>(face_index[3U*i]);
        const auto neighbour = face_index[3U*i+1U] == invalid ? cfd::fvm::invalid_cell : static_cast<std::size_t>(face_index[3U*i+1U]);
        const auto patch = face_index[3U*i+2U] == invalid ? cfd::fvm::invalid_patch : static_cast<std::size_t>(face_index[3U*i+2U]);
        faces.push_back({owner, neighbour,
                         {face_geometry[6U*i], face_geometry[6U*i+1U], face_geometry[6U*i+2U]},
                         {face_geometry[6U*i+3U], face_geometry[6U*i+4U], face_geometry[6U*i+5U]}, patch});
    }

    const auto patch_names = read_string_table(file.get(), "/Mesh/PatchNames");
    std::vector<cfd::fvm::BoundaryPatch> patches;
    patches.reserve(patch_names.size());
    for (const auto& name : patch_names) patches.push_back({name});
    cfd::fvm::PolyMesh mesh(std::move(cells), std::move(faces), std::move(patches));

    const hsize_t one[1]{1U};
    const auto time_value = read_double_dataset(file.get(), "/Meta/Time", one);
    const auto step_value = read_u64_dataset(file.get(), "/Meta/Step", one);
    require_h5(std::isfinite(time_value.front()) && time_value.front() >= 0.0, "invalid FVM checkpoint time");
    const auto field_names = read_string_table(file.get(), "/Meta/FieldNames");
    std::vector<FvmScalarField> fields;
    fields.reserve(field_names.size());
    const hsize_t field_dims[1]{static_cast<hsize_t>(mesh.cell_count())};
    for (const auto& name : field_names) {
        validate_field_name(name);
        fields.push_back({name, read_double_dataset(file.get(), "/Fields/" + name, field_dims)});
    }
    return {std::move(mesh), time_value.front(), static_cast<std::size_t>(step_value.front()), std::move(fields)};
#endif
}

} // namespace cfd::io

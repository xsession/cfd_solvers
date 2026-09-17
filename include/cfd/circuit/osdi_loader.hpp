#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

namespace cfd::circuit {

// Minimal dynamic-library seam for Open Source Device Interface (OSDI)
// modules produced by OpenVAF. This loader intentionally stops at ABI
// discovery/validation: evaluation of OsdiDescriptor instances remains a
// separate adapter layer so the core simulator does not depend on OSDI's
// evolving descriptor layout.
class OsdiLibrary {
public:
    explicit OsdiLibrary(const std::string& path);
    ~OsdiLibrary();

    OsdiLibrary(const OsdiLibrary&) = delete;
    OsdiLibrary& operator=(const OsdiLibrary&) = delete;
    OsdiLibrary(OsdiLibrary&& other) noexcept;
    OsdiLibrary& operator=(OsdiLibrary&& other) noexcept;

    [[nodiscard]] std::uint32_t version_major() const noexcept { return version_major_; }
    [[nodiscard]] std::uint32_t version_minor() const noexcept { return version_minor_; }
    [[nodiscard]] std::uint32_t descriptor_count() const noexcept { return descriptor_count_; }
    [[nodiscard]] const void* raw_descriptors() const noexcept { return descriptors_; }
    [[nodiscard]] const std::string& path() const noexcept { return path_; }

private:
    void* handle_{};
    std::uint32_t version_major_{};
    std::uint32_t version_minor_{};
    std::uint32_t descriptor_count_{};
    const void* descriptors_{};
    std::string path_;

    void close() noexcept;
};

} // namespace cfd::circuit

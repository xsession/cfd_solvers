#pragma once
#include <string_view>
#include <vector>
namespace cfd::core {
enum class RuntimeBackend { serial, openmp, mpi, sycl };
struct BackendInfo { RuntimeBackend backend{}; const char* name{}; bool available{}; };
[[nodiscard]] std::vector<BackendInfo> compiled_backends();
[[nodiscard]] RuntimeBackend parse_backend(std::string_view name);
[[nodiscard]] RuntimeBackend select_backend(std::string_view request="auto");
[[nodiscard]] const char* backend_name(RuntimeBackend backend) noexcept;
[[nodiscard]] bool backend_available(RuntimeBackend backend) noexcept;
} // namespace cfd::core

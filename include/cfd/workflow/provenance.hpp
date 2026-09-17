#pragma once
#include <filesystem>
#include <map>
#include <string>
namespace cfd::workflow {
struct ResultProvenance {std::string solver;std::string version;std::string git_commit;std::string timestamp_utc;std::map<std::string,std::string> metadata;};
[[nodiscard]] ResultProvenance make_provenance(std::string solver,std::string version,std::string git_commit="unknown");
[[nodiscard]] std::string provenance_json(const ResultProvenance& p);void write_provenance(const std::filesystem::path& path,const ResultProvenance& p);
} // namespace cfd::workflow

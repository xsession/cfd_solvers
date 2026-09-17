#include "cfd/core/autotune.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace cfd::core {
namespace {
std::string clean_field(std::string value) {
    for (char& ch : value) {
        if (ch == '\t' || ch == '\n' || ch == '\r') ch = ' ';
    }
    return value;
}

std::vector<std::string> split_tsv(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0U;
    while (true) {
        const auto pos = line.find('\t', begin);
        if (pos == std::string::npos) {
            fields.push_back(line.substr(begin));
            break;
        }
        fields.push_back(line.substr(begin, pos - begin));
        begin = pos + 1U;
    }
    return fields;
}

std::size_t parse_size(const std::string& text, const char* what) {
    std::size_t used = 0U;
    const auto value = std::stoull(text, &used);
    if (used != text.size() || value == 0U || value > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error(std::string("invalid tuning ") + what);
    }
    return static_cast<std::size_t>(value);
}
} // namespace

std::string device_tuning_key_string(const DeviceTuningKey& key) {
    return clean_field(key.vendor) + "|" + clean_field(key.name) + "|" + clean_field(key.driver);
}

KernelTuningDatabase KernelTuningDatabase::load(const std::filesystem::path& path) {
    KernelTuningDatabase database;
    std::ifstream input(path);
    if (!input) return database;
    std::string line;
    std::size_t line_number = 0U;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') continue;
        const auto fields = split_tsv(line);
        if (fields.size() != 8U) {
            throw std::runtime_error("invalid tuning database row at line " + std::to_string(line_number));
        }
        KernelTuningRecord record;
        record.device = {fields[0], fields[1], fields[2]};
        record.kernel_group = fields[3];
        record.candidate.work_group_size = parse_size(fields[4], "work-group size");
        record.candidate.vector_width = parse_size(fields[5], "vector width");
        record.candidate.fused_stages = parse_size(fields[6], "fused stage count");
        std::size_t used = 0U;
        record.score_seconds = std::stod(fields[7], &used);
        if (used != fields[7].size() || !(record.score_seconds > 0.0) || !std::isfinite(record.score_seconds)) {
            throw std::runtime_error("invalid tuning score at line " + std::to_string(line_number));
        }
        database.upsert(std::move(record));
    }
    return database;
}

void KernelTuningDatabase::save(const std::filesystem::path& path) const {
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create tuning database: " + path.string());
    output << "# vendor\tdevice\tdriver\tkernel_group\twork_group\tvector_width\tfused_stages\tscore_seconds\n";
    output.precision(17);
    for (const auto& record : records_) {
        output << clean_field(record.device.vendor) << '\t'
               << clean_field(record.device.name) << '\t'
               << clean_field(record.device.driver) << '\t'
               << clean_field(record.kernel_group) << '\t'
               << record.candidate.work_group_size << '\t'
               << record.candidate.vector_width << '\t'
               << record.candidate.fused_stages << '\t'
               << record.score_seconds << '\n';
    }
}

void KernelTuningDatabase::upsert(KernelTuningRecord record) {
    if (record.kernel_group.empty() || record.candidate.work_group_size == 0U ||
        record.candidate.vector_width == 0U || record.candidate.fused_stages == 0U ||
        !(record.score_seconds > 0.0) || !std::isfinite(record.score_seconds)) {
        throw std::invalid_argument("invalid kernel tuning record");
    }
    const auto match = [&](const KernelTuningRecord& existing) {
        return existing.device == record.device && existing.kernel_group == record.kernel_group;
    };
    const auto it = std::find_if(records_.begin(), records_.end(), match);
    if (it == records_.end()) {
        records_.push_back(std::move(record));
    } else if (record.score_seconds <= it->score_seconds) {
        *it = std::move(record);
    }
}

std::optional<KernelTuningRecord> KernelTuningDatabase::lookup(const DeviceTuningKey& device,
                                                                std::string_view kernel_group) const {
    const auto it = std::find_if(records_.begin(), records_.end(), [&](const KernelTuningRecord& record) {
        return record.device == device && record.kernel_group == kernel_group;
    });
    if (it == records_.end()) return std::nullopt;
    return *it;
}

#if defined(CFD_HAS_SYCL)
DeviceTuningKey sycl_device_tuning_key(const sycl::device& device) {
    return {
        device.get_info<sycl::info::device::vendor>(),
        device.get_info<sycl::info::device::name>(),
        device.get_info<sycl::info::device::driver_version>()
    };
}
#endif

} // namespace cfd::core

#include "cfd/chemistry/database_io.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace cfd::chemistry {
namespace {

[[nodiscard]] std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r");
    if (first == std::string::npos)
        return {};
    const auto last = value.find_last_not_of(" \t\r");
    return value.substr(first, last - first + 1U);
}

[[nodiscard]] std::string lower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

[[nodiscard]] double numeric_after(const std::string& line, std::string_view key, double fallback) {
    const auto found = lower(line).find(lower(std::string(key)));
    if (found == std::string::npos)
        return fallback;
    const auto colon = line.find(':', found);
    const auto begin = colon == std::string::npos ? found + key.size() : colon + 1U;
    std::istringstream in(line.substr(begin));
    double value = fallback;
    in >> value;
    return in ? value : fallback;
}

[[nodiscard]] std::size_t find_species(const std::unordered_map<std::string, std::size_t>& indices, std::string name) {
    name = trim(std::move(name));
    while (!name.empty() && (name.front() == '+' || name.front() == '-'))
        name.erase(name.begin());
    const auto found = indices.find(name);
    if (found == indices.end())
        throw std::invalid_argument("Cantera reaction references unknown species: " + name);
    return found->second;
}

[[nodiscard]] std::vector<StoichiometricTerm> parse_side(std::string side,
                                                         const std::unordered_map<std::string, std::size_t>& indices) {
    std::vector<StoichiometricTerm> terms;
    std::size_t begin = 0U;
    while (begin < side.size()) {
        const auto plus = side.find('+', begin);
        std::string term = trim(side.substr(begin, plus == std::string::npos ? std::string::npos : plus - begin));
        if (!term.empty()) {
            std::istringstream in(term);
            double coefficient = 1.0;
            std::string name;
            if (!(in >> name))
                throw std::invalid_argument("empty Cantera stoichiometric term");
            try {
                std::size_t used = 0U;
                const double parsed = std::stod(name, &used);
                if (used == name.size()) {
                    coefficient = parsed;
                    if (!(in >> name))
                        throw std::invalid_argument("missing Cantera species name");
                } else if (used > 0U) {
                    coefficient = parsed;
                    name = name.substr(used);
                }
            } catch (const std::exception&) {
            }
            if (!(coefficient > 0.0) || !std::isfinite(coefficient))
                throw std::invalid_argument("invalid Cantera stoichiometric coefficient");
            terms.push_back({find_species(indices, std::move(name)), coefficient});
        }
        if (plus == std::string::npos)
            break;
        begin = plus + 1U;
    }
    return terms;
}

void add_equation(ReactionNetwork& network, const std::unordered_map<std::string, std::size_t>& indices,
                  const std::string& equation, double pre_exponential, double temperature_exponent,
                  double activation_energy) {
    const auto reversible = equation.find("<=>");
    const auto irreversible = equation.find("=>");
    const auto arrow = reversible != std::string::npos ? reversible : irreversible;
    if (arrow == std::string::npos)
        throw std::invalid_argument("Cantera reaction is missing an arrow");
    const std::size_t arrow_width = reversible != std::string::npos ? 3U : 2U;
    ElementaryReaction reaction;
    reaction.reactants = parse_side(equation.substr(0U, arrow), indices);
    reaction.products = parse_side(equation.substr(arrow + arrow_width), indices);
    reaction.forward = {pre_exponential, temperature_exponent, activation_energy};
    network.add_reaction(std::move(reaction));
}

} // namespace

std::vector<Species> load_phreeqc_master_species(std::string_view database_text) {
    std::istringstream input{std::string(database_text)};
    std::vector<Species> species;
    std::string line;
    bool in_section = false;
    while (std::getline(input, line)) {
        const std::string clean = trim(line.substr(0U, line.find('#')));
        if (clean.empty())
            continue;
        if (lower(clean) == "solution_master_species") {
            in_section = true;
            continue;
        }
        if (in_section && std::isupper(static_cast<unsigned char>(clean.front())) &&
            clean.find_first_of(" \t") == std::string::npos) {
            in_section = false;
            continue;
        }
        if (!in_section)
            continue;
        std::istringstream row(clean);
        std::string element, name;
        double molar_mass = 0.0;
        double charge = 0.0;
        if (!(row >> element >> name >> molar_mass >> charge))
            continue;
        if (!(molar_mass > 0.0) || !std::isfinite(molar_mass) || !std::isfinite(charge))
            throw std::invalid_argument("invalid PHREEQC master species row");
        if (molar_mass > 1.0)
            molar_mass *= 1.0e-3;
        species.push_back({std::move(name), molar_mass, static_cast<int>(std::lround(charge))});
    }
    if (species.empty())
        throw std::invalid_argument("PHREEQC database has no SOLUTION_MASTER_SPECIES rows");
    return species;
}

ReactionNetwork load_cantera_yaml(std::string_view mechanism_text) {
    std::istringstream input{std::string(mechanism_text)};
    std::vector<Species> species;
    std::unordered_map<std::string, std::size_t> indices;
    std::string line;
    bool in_species = false;
    bool in_reactions = false;
    std::string equation;
    double pre_exponential = 0.0;
    double temperature_exponent = 0.0;
    double activation_energy = 0.0;
    ReactionNetwork network;
    const auto flush = [&]() {
        if (!equation.empty()) {
            add_equation(network, indices, equation, pre_exponential, temperature_exponent, activation_energy);
            equation.clear();
            pre_exponential = 0.0;
            temperature_exponent = 0.0;
            activation_energy = 0.0;
        }
    };
    while (std::getline(input, line)) {
        const std::string clean = trim(line.substr(0U, line.find('#')));
        if (clean.empty())
            continue;
        const std::string folded = lower(clean);
        if (folded == "species:") {
            in_species = true;
            in_reactions = false;
            continue;
        }
        if (folded == "reactions:") {
            in_species = false;
            in_reactions = true;
            continue;
        }
        if (in_species && folded.rfind("- name:", 0U) == 0U) {
            std::string name = trim(clean.substr(clean.find(':') + 1U));
            if (name.empty())
                throw std::invalid_argument("Cantera species name is empty");
            if (indices.contains(name))
                throw std::invalid_argument("duplicate Cantera species: " + name);
            indices[name] = species.size();
            species.push_back({name, 1.0, 0});
            (void)network.add_species(species.back());
            continue;
        }
        if (in_reactions && folded.rfind("- equation:", 0U) == 0U) {
            flush();
            equation = trim(clean.substr(clean.find(':') + 1U));
            continue;
        }
        if (in_reactions && !equation.empty()) {
            pre_exponential = numeric_after(clean, "a:", pre_exponential);
            temperature_exponent = numeric_after(clean, "b:", temperature_exponent);
            activation_energy = numeric_after(clean, "ea:", activation_energy);
        }
    }
    flush();
    if (species.empty())
        throw std::invalid_argument("Cantera YAML has no species section");
    return network;
}

ReaktoroDynamicAdapter::ReaktoroDynamicAdapter(Evaluator evaluator) : evaluator_(std::move(evaluator)) {
    if (!evaluator_)
        throw std::invalid_argument("Reaktoro adapter evaluator must not be empty");
}

std::vector<double> ReaktoroDynamicAdapter::equilibrate(std::span<const double> composition, double temperature,
                                                        double pressure) const {
    if (!(temperature > 0.0) || !(pressure > 0.0) || !std::isfinite(temperature) || !std::isfinite(pressure)) {
        throw std::invalid_argument("invalid Reaktoro adapter state");
    }
    auto result = evaluator_(composition, temperature, pressure);
    if (result.size() != composition.size())
        throw std::runtime_error("Reaktoro adapter size mismatch");
    for (double value : result)
        if (!(value >= 0.0) || !std::isfinite(value))
            throw std::runtime_error("Reaktoro adapter returned invalid composition");
    return result;
}

} // namespace cfd::chemistry

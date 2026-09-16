#pragma once

#include "cfd/core/conjugate_gradient.hpp"
#include "cfd/fvm/poly_mesh.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace cfd::electrochemistry {

enum class SpeciesBoundaryType {
    noFlux,
    fixedConcentration,
    fixedMolarFlux,
    butlerVolmerFlux,
};

enum class PotentialBoundaryType {
    insulating,
    fixedPotential,
};

struct PotentialBoundaryCondition {
    PotentialBoundaryType type{PotentialBoundaryType::insulating};
    double value{};
};

struct SpeciesBoundaryCondition {
    SpeciesBoundaryType type{SpeciesBoundaryType::noFlux};
    double value{}; // concentration [mol/m3] or outward molar flux [mol/m2/s]

    // Butler-Volmer data. Positive anodic current produces `stoichiometric`
    // moles of this species per reaction event into the electrolyte.
    double electrode_potential{};
    double equilibrium_potential{};
    double exchange_current_density{};
    double alpha_anodic{0.5};
    double alpha_cathodic{0.5};
    int electrons{1};
    double stoichiometric{1.0};
};

struct PolyTransportSpecies {
    std::string name;
    int charge{};
    double diffusivity{};
    std::vector<double> concentration;
    std::vector<SpeciesBoundaryCondition> boundary;
};

enum class ElectromigrationFluxScheme {
    centered,
    scharfetterGummel,
};

struct NernstPlanckPolyConfig {
    double temperature{298.15};
    double dt{1.0e-5};
    double relative_permittivity{78.5};
    ElectromigrationFluxScheme electromigration_scheme{ElectromigrationFluxScheme::centered};
};

// Explicit conservative finite-volume Nernst-Planck transport on PolyMesh.
// Potential and advective volume flux are supplied by the surrounding
// electrostatic/flow solver, allowing this transport kernel to be reused by
// electroneutral, Poisson-Nernst-Planck and coupled CFD formulations.
class NernstPlanckPolyMesh {
public:
    explicit NernstPlanckPolyMesh(cfd::fvm::PolyMesh mesh, NernstPlanckPolyConfig config = {});

    std::size_t add_species(std::string name,
                            int charge,
                            double diffusivity,
                            double initial_concentration);
    std::size_t add_species(std::string name,
                            int charge,
                            double diffusivity,
                            const std::function<double(cfd::fvm::Vec3)>& initial_concentration);
    void set_species_boundary(std::size_t species,
                              std::string_view patch,
                              SpeciesBoundaryType type,
                              double value = 0.0);
    void set_species_butler_volmer_boundary(std::size_t species,
                                            std::string_view patch,
                                            double electrode_potential,
                                            double equilibrium_potential,
                                            double exchange_current_density,
                                            int electrons,
                                            double stoichiometric = 1.0,
                                            double alpha_anodic = 0.5,
                                            double alpha_cathodic = 0.5);
    void set_potential(double value);
    void set_potential(const std::function<double(cfd::fvm::Vec3)>& value);
    void set_potential_boundary(std::string_view patch,
                                PotentialBoundaryType type,
                                double value = 0.0);
    void set_face_flux(std::vector<double> volumetric_flux);

    void step(std::size_t steps = 1U);
    void step_poisson_nernst_planck(std::size_t steps = 1U,
                                    std::size_t poisson_max_iterations = 1000U,
                                    double poisson_relative_tolerance = 1.0e-10);

    // Dilute-solution electroneutral potential: enforce div(i)=0 with
    // i = -F sum(z D grad(c)) - kappa grad(phi). All boundaries are
    // insulating/current-free in this baseline; the potential gauge is mean-zero.
    cfd::core::ConjugateGradientResult solve_electroneutral_potential(
        std::size_t max_iterations = 1000U, double relative_tolerance = 1.0e-10);

    // Poisson-Nernst-Planck electrostatics:
    //   -div(epsilon grad(phi)) = F sum(z_k c_k).
    // Fixed-potential and insulating patch conditions are supported. A pure
    // insulating problem requires net-zero charge and uses a mean-zero gauge.
    cfd::core::ConjugateGradientResult solve_poisson_potential(
        std::size_t max_iterations = 1000U, double relative_tolerance = 1.0e-10);
    [[nodiscard]] std::vector<double> current_flux_faces() const;
    [[nodiscard]] double total_charge() const;

    [[nodiscard]] const cfd::fvm::PolyMesh& mesh() const noexcept { return mesh_; }
    [[nodiscard]] const std::vector<PolyTransportSpecies>& species() const noexcept { return species_; }
    [[nodiscard]] const std::vector<double>& potential() const noexcept { return potential_; }
    [[nodiscard]] double time() const noexcept { return time_; }
    [[nodiscard]] double total_amount(std::size_t species) const;
    [[nodiscard]] double minimum_concentration(std::size_t species) const;

private:
    cfd::fvm::PolyMesh mesh_;
    NernstPlanckPolyConfig config_;
    std::vector<PolyTransportSpecies> species_;
    std::vector<double> potential_;
    std::vector<PotentialBoundaryCondition> potential_boundary_;
    std::vector<double> face_flux_;
    cfd::core::ConjugateGradientWorkspace potential_workspace_;
    cfd::core::ConjugateGradientResult potential_result_{};
    double time_{};

    void step_once();
};

} // namespace cfd::electrochemistry

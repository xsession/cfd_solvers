#pragma once
#include <array>
namespace cfd::fem {
struct MaterialPoint1D {double stress{},tangent{},plastic_strain{},equivalent_plastic_strain{};};
class UniaxialMaterialLaw {public:virtual ~UniaxialMaterialLaw()=default;virtual MaterialPoint1D update(double total_strain)=0;};
class LinearElastic1D final:public UniaxialMaterialLaw {public:explicit LinearElastic1D(double youngs);MaterialPoint1D update(double total_strain)override;private:double E_;};
class BilinearPlasticity1D final:public UniaxialMaterialLaw {public:BilinearPlasticity1D(double youngs,double yield_stress,double hardening_modulus);MaterialPoint1D update(double total_strain)override;private:double E_,yield_,H_,plastic_{},alpha_{};};
struct NeoHookeanResult {double energy_density{};std::array<double,3> nominal_stress{};};
[[nodiscard]] NeoHookeanResult compressible_neo_hookean(std::array<double,3> principal_stretch,double shear_modulus,double lame_lambda);
[[nodiscard]] double penalty_contact_pressure(double signed_gap,double penalty_stiffness);
}

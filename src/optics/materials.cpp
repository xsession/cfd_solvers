#include "cfd/solvers/optics/materials.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::optics {

double SellmeierMaterial::refractive_index_nm(double wavelength_nm) const{
    if(!(wavelength_nm>0.0)) throw std::invalid_argument("optical wavelength must be positive");
    const double lambda_um=1.0e-3*wavelength_nm;
    const double l2=lambda_um*lambda_um;
    double n2=1.0;
    for(std::size_t i=0;i<3U;++i){
        const double denom=l2-c_um2[i];
        if(std::abs(denom)<1.0e-14) throw std::domain_error("Sellmeier wavelength at material resonance");
        n2+=b[i]*l2/denom;
    }
    if(!(n2>0.0)) throw std::domain_error("Sellmeier model produced non-positive n^2");
    return std::sqrt(n2);
}

MaterialCatalog::MaterialCatalog(){
    // Common public Sellmeier forms. C coefficients are in micrometre^2.
    add({"N-BK7",{1.03961212,0.231792344,1.01046945},{0.00600069867,0.0200179144,103.560653}});
    add({"FusedSilica",{0.6961663,0.4079426,0.8974794},{0.0684043*0.0684043,0.1162414*0.1162414,9.896161*9.896161}});
}
void MaterialCatalog::add(SellmeierMaterial material){
    if(material.name.empty()) throw std::invalid_argument("optical material name must not be empty");
    if(std::any_of(materials_.begin(),materials_.end(),[&](const auto&m){return m.name==material.name;})) throw std::invalid_argument("duplicate optical material name");
    materials_.push_back(std::move(material));
}
const SellmeierMaterial& MaterialCatalog::at(std::string_view name) const{
    const auto it=std::find_if(materials_.begin(),materials_.end(),[&](const auto&m){return m.name==name;});
    if(it==materials_.end()) throw std::out_of_range("unknown optical material");
    return *it;
}


const SellmeierMaterial& MaterialCatalog::nearest_index(double wavelength_nm,double target_index) const{
    if(materials_.empty()||!(target_index>0.0)||!std::isfinite(target_index))throw std::invalid_argument("invalid material search");
    const SellmeierMaterial* best=&materials_.front();double best_error=std::abs(best->refractive_index_nm(wavelength_nm)-target_index);
    for(const auto& material:materials_){const double error=std::abs(material.refractive_index_nm(wavelength_nm)-target_index);if(error<best_error){best=&material;best_error=error;}}
    return *best;
}
} // namespace cfd::optics

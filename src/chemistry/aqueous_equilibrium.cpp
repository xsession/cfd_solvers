#include "cfd/chemistry/aqueous_equilibrium.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace cfd::chemistry {
namespace {
std::vector<double> distribution(const AcidFamily& family,double log_h) {
    std::vector<double> weights(family.dissociation_constants.size()+1U,0.0);
    for(std::size_t i=1;i<weights.size();++i)
        weights[i]=weights[i-1]+std::log(family.dissociation_constants[i-1])-log_h;
    const double maximum=*std::max_element(weights.begin(),weights.end());
    double sum=0.0;
    for(double& w:weights){w=std::exp(w-maximum);sum+=w;}
    for(double& w:weights)w*=family.total_mol_per_litre/sum;
    return weights;
}
}

AqueousEquilibriumResult equilibrate_acids(std::span<const AcidFamily> families,double strong_charge,double kw) {
    if(!std::isfinite(strong_charge)||!(kw>0.0)||!std::isfinite(kw))
        throw std::invalid_argument("invalid aqueous equilibrium controls");
    for(const auto& family:families){
        if(!(family.total_mol_per_litre>=0.0)||!std::isfinite(family.total_mol_per_litre))
            throw std::invalid_argument("invalid analytical acid concentration");
        for(double ka:family.dissociation_constants)if(!(ka>0.0)||!std::isfinite(ka))
            throw std::invalid_argument("invalid acid dissociation constant");
    }
    const auto charge=[&](double log_h){
        double q=std::exp(log_h)-kw/std::exp(log_h)+strong_charge;
        for(const auto& family:families){
            const auto concentration=distribution(family,log_h);
            for(std::size_t j=0;j<concentration.size();++j)
                q+=(static_cast<double>(family.fully_protonated_charge)-static_cast<double>(j))*concentration[j];
        }
        return q;
    };
    double low=-690.0,high=690.0;
    if(!(charge(low)<0.0&&charge(high)>0.0)) throw std::runtime_error("aqueous charge balance cannot be bracketed");
    for(std::size_t i=0;i<200U;++i){
        const double mid=0.5*(low+high);
        if(charge(mid)>0.0)high=mid;else low=mid;
        if(high-low<1.0e-13)break;
    }
    const double log_h=0.5*(low+high),h=std::exp(log_h);
    AqueousEquilibriumResult result{-log_h/std::log(10.0),h,kw/h,charge(log_h),{}};
    for(const auto& family:families)result.family_concentrations.push_back(distribution(family,log_h));
    return result;
}
} // namespace cfd::chemistry

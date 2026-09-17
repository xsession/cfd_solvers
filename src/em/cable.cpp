#include "cfd/em/cable.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace cfd::em {
namespace {
using Matrix=cfd::rf::ComplexMatrix;
Matrix add(Matrix a,const Matrix& b,double scale=1.0){if(a.size()!=b.size())throw std::invalid_argument("cable matrix size mismatch");for(std::size_t r=0;r<a.size();++r)for(std::size_t c=0;c<a.size();++c)a(r,c)+=scale*b(r,c);return a;}
std::vector<Complex> matvec(const Matrix& a,std::span<const Complex> x){if(a.size()!=x.size())throw std::invalid_argument("cable matrix-vector size mismatch");std::vector<Complex> y(x.size());for(std::size_t r=0;r<a.size();++r)for(std::size_t c=0;c<a.size();++c)y[r]+=a(r,c)*x[c];return y;}
Matrix block(const Matrix& source,std::size_t row0,std::size_t col0,std::size_t n){Matrix out(n);for(std::size_t r=0;r<n;++r)for(std::size_t c=0;c<n;++c)out(r,c)=source(row0+r,col0+c);return out;}
Matrix subtract_product(const Matrix& a,const Matrix& b,const Matrix& c){return add(a,cfd::rf::multiply(b,c),-1.0);}
}

void MulticonductorRlcg::validate() const {if(conductors==0U)throw std::invalid_argument("cable RLCG requires conductors");const std::size_t n=conductors*conductors;if(resistance.size()!=n||inductance.size()!=n||conductance.size()!=n||capacitance.size()!=n)throw std::invalid_argument("cable RLCG matrix size mismatch");for(double value:resistance)if(!std::isfinite(value))throw std::invalid_argument("non-finite cable resistance");for(double value:inductance)if(!std::isfinite(value))throw std::invalid_argument("non-finite cable inductance");for(double value:conductance)if(!std::isfinite(value))throw std::invalid_argument("non-finite cable conductance");for(double value:capacitance)if(!std::isfinite(value))throw std::invalid_argument("non-finite cable capacitance");}

CableFrequencyResult solve_multiconductor_cable(double length,std::size_t segments,double frequency,const MulticonductorRlcg& rlcg,std::span<const Complex> input_voltage,const Matrix& load){
    rlcg.validate();const std::size_t n=rlcg.conductors;if(!(length>0.0)||segments==0U||!(frequency>=0.0)||input_voltage.size()!=n||load.size()!=n)throw std::invalid_argument("invalid multiconductor cable controls");const double omega=2.0*std::numbers::pi*frequency,dz=length/static_cast<double>(segments);Matrix system(2U*n);for(std::size_t r=0;r<n;++r)for(std::size_t c=0;c<n;++c){const std::size_t k=r*n+c;const Complex z{rlcg.resistance[k],omega*rlcg.inductance[k]},y{rlcg.conductance[k],omega*rlcg.capacitance[k]};system(r,n+c)=-z;system(n+r,c)=-y;}
    Matrix lhs=Matrix::identity(2U*n),rhs=Matrix::identity(2U*n);for(std::size_t r=0;r<2U*n;++r)for(std::size_t c=0;c<2U*n;++c){lhs(r,c)-=0.5*dz*system(r,c);rhs(r,c)+=0.5*dz*system(r,c);}const Matrix step=cfd::rf::multiply(cfd::rf::inverse(lhs),rhs);Matrix transfer=Matrix::identity(2U*n);for(std::size_t k=0;k<segments;++k)transfer=cfd::rf::multiply(step,transfer);
    const Matrix a=block(transfer,0U,0U,n),b=block(transfer,0U,n,n),c=block(transfer,n,0U,n),d=block(transfer,n,n,n);const Matrix left=subtract_product(b,load,d);const Matrix right=add(cfd::rf::multiply(load,c),a,-1.0);const auto right_v=matvec(right,input_voltage);const auto input_current=matvec(cfd::rf::inverse(left),right_v);std::vector<Complex> state(2U*n);for(std::size_t i=0;i<n;++i){state[i]=input_voltage[i];state[n+i]=input_current[i];}const auto load_state=matvec(transfer,state);CableFrequencyResult result;result.frequency_hz=frequency;result.input_voltage_v.assign(input_voltage.begin(),input_voltage.end());result.input_current_a=input_current;result.load_voltage_v.assign(load_state.begin(),load_state.begin()+static_cast<std::ptrdiff_t>(n));result.load_current_a.assign(load_state.begin()+static_cast<std::ptrdiff_t>(n),load_state.end());result.transfer_matrix=transfer;return result;
}

std::vector<CableFrequencyResult> sweep_multiconductor_cable(double length,std::size_t segments,std::span<const double> frequencies,const RlcgSampler& sampler,std::span<const Complex> input_voltage,const Matrix& load){if(frequencies.empty()||!sampler)throw std::invalid_argument("invalid cable frequency sweep");std::vector<CableFrequencyResult> out;out.reserve(frequencies.size());for(double frequency:frequencies)out.push_back(solve_multiconductor_cable(length,segments,frequency,sampler(frequency),input_voltage,load));return out;}

Complex cable_shield_transfer_impedance_ohm_per_m(double frequency,const CableShieldTransferModel& model){if(!(frequency>=0.0)||!(model.dc_transfer_resistance_ohm_per_m>=0.0)||!(model.transfer_inductance_h_per_m>=0.0)||!(model.skin_corner_hz>0.0))throw std::invalid_argument("invalid cable shield transfer model");const Complex skin=std::sqrt(Complex{1.0,frequency/model.skin_corner_hz});return model.dc_transfer_resistance_ohm_per_m*skin+Complex{0.0,2.0*std::numbers::pi*frequency*model.transfer_inductance_h_per_m};}
Complex cable_shield_induced_voltage_v(double frequency,Complex current,double length,const CableShieldTransferModel& model){if(!(length>=0.0))throw std::invalid_argument("invalid cable exposed length");return current*cable_shield_transfer_impedance_ohm_per_m(frequency,model)*length;}
Complex field_to_cable_open_circuit_voltage_v(Complex field,double effective_height,Complex coupling){if(!(effective_height>=0.0)||!std::isfinite(field.real())||!std::isfinite(field.imag())||!std::isfinite(coupling.real())||!std::isfinite(coupling.imag()))throw std::invalid_argument("invalid field-to-cable coupling controls");return field*effective_height*coupling;}

} // namespace cfd::em

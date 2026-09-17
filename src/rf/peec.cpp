#include "cfd/rf/peec.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace cfd::rf {
namespace {
constexpr double mu0 = 4.0e-7 * std::numbers::pi;
constexpr double eps0 = 8.8541878128e-12;

Point3 sub(Point3 a, Point3 b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
Point3 add(Point3 a, Point3 b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
Point3 scale(Point3 a, double s) { return {a.x*s, a.y*s, a.z*s}; }
double dot(Point3 a, Point3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
double norm(Point3 a) { return std::sqrt(dot(a,a)); }

struct Rule { std::span<const double> x, w; };
Rule gauss_rule(std::size_t order) {
    static constexpr std::array<double,2> x2{-0.57735026918962576451, 0.57735026918962576451};
    static constexpr std::array<double,2> w2{1.0, 1.0};
    static constexpr std::array<double,4> x4{-0.86113631159405257522,-0.33998104358485626480,
                                              0.33998104358485626480, 0.86113631159405257522};
    static constexpr std::array<double,4> w4{0.34785484513745385737,0.65214515486254614263,
                                              0.65214515486254614263,0.34785484513745385737};
    static constexpr std::array<double,6> x6{-0.93246951420315202781,-0.66120938646626451366,
                                             -0.23861918608319690863, 0.23861918608319690863,
                                              0.66120938646626451366, 0.93246951420315202781};
    static constexpr std::array<double,6> w6{0.17132449237917034504,0.36076157304813860757,
                                              0.46791393457269104739,0.46791393457269104739,
                                              0.36076157304813860757,0.17132449237917034504};
    if (order <= 2U) return {x2,w2};
    if (order <= 4U) return {x4,w4};
    return {x6,w6};
}

template<class T>
std::vector<T> solve_dense(std::vector<T> matrix, std::vector<T> rhs, std::size_t n) {
    for (std::size_t column=0; column<n; ++column) {
        std::size_t pivot=column;
        double best=std::abs(matrix[column*n+column]);
        for (std::size_t row=column+1; row<n; ++row) {
            if (const double candidate=std::abs(matrix[row*n+column]); candidate>best) {
                best=candidate; pivot=row;
            }
        }
        if (best<1.0e-24) throw std::runtime_error("singular PEEC impedance matrix");
        if (pivot!=column) {
            for (std::size_t j=column; j<n; ++j) std::swap(matrix[column*n+j],matrix[pivot*n+j]);
            std::swap(rhs[column],rhs[pivot]);
        }
        const T diagonal=matrix[column*n+column];
        for (std::size_t j=column; j<n; ++j) matrix[column*n+j]/=diagonal;
        rhs[column]/=diagonal;
        for (std::size_t row=0; row<n; ++row) {
            if (row==column) continue;
            const T factor=matrix[row*n+column];
            for (std::size_t j=column; j<n; ++j) matrix[row*n+j]-=factor*matrix[column*n+j];
            rhs[row]-=factor*rhs[column];
        }
    }
    return rhs;
}

std::vector<double> inverse_dense_real(std::vector<double> matrix, std::size_t n) {
    std::vector<double> inverse(n*n,0.0);
    for(std::size_t i=0;i<n;++i) inverse[i*n+i]=1.0;
    for(std::size_t column=0;column<n;++column){
        std::size_t pivot=column;double best=std::abs(matrix[column*n+column]);
        for(std::size_t row=column+1U;row<n;++row){
            const double candidate=std::abs(matrix[row*n+column]);
            if(candidate>best){best=candidate;pivot=row;}
        }
        if(best<1.0e-30)throw std::runtime_error("singular PEEC potential matrix");
        if(pivot!=column){
            for(std::size_t j=0;j<n;++j){
                std::swap(matrix[column*n+j],matrix[pivot*n+j]);
                std::swap(inverse[column*n+j],inverse[pivot*n+j]);
            }
        }
        const double diagonal=matrix[column*n+column];
        for(std::size_t j=0;j<n;++j){matrix[column*n+j]/=diagonal;inverse[column*n+j]/=diagonal;}
        for(std::size_t row=0;row<n;++row){
            if(row==column)continue;
            const double factor=matrix[row*n+column];
            if(factor==0.0)continue;
            for(std::size_t j=0;j<n;++j){
                matrix[row*n+j]-=factor*matrix[column*n+j];
                inverse[row*n+j]-=factor*inverse[column*n+j];
            }
        }
    }
    return inverse;
}

double potential_coefficient(const FilamentSegment& first,const FilamentSegment& second,
                             double relative_permittivity,std::size_t order,bool self) {
    const Point3 d1=sub(first.end,first.start),d2=sub(second.end,second.start);
    const auto rule=gauss_rule(order);double average_kernel=0.0;
    for(std::size_t i=0;i<rule.x.size();++i){
        const double u=0.5*(rule.x[i]+1.0);const Point3 p1=add(first.start,scale(d1,u));
        for(std::size_t j=0;j<rule.x.size();++j){
            const double v=0.5*(rule.x[j]+1.0);const Point3 p2=add(second.start,scale(d2,v));
            const Point3 delta=sub(p1,p2);
            const double smoothing=self?first.radius_m*first.radius_m:0.0;
            const double distance=std::sqrt(dot(delta,delta)+smoothing);
            if(!(distance>0.0))throw std::runtime_error("singular PEEC coefficient-of-potential integral");
            average_kernel+=0.25*rule.w[i]*rule.w[j]/distance;
        }
    }
    return average_kernel/(4.0*std::numbers::pi*eps0*relative_permittivity);
}

void validate(const FilamentSegment& segment) {
    const double length=filament_length(segment);
    if (!(length>0.0) || !(segment.radius_m>0.0) || !(segment.conductivity_s_per_m>0.0)
        || !std::isfinite(segment.radius_m) || !std::isfinite(segment.conductivity_s_per_m)) {
        throw std::invalid_argument("invalid PEEC filament segment");
    }
    if (segment.radius_m>=0.25*length) {
        throw std::invalid_argument("PEEC filament requires a slender conductor");
    }
}
} // namespace

double filament_length(const FilamentSegment& segment) {
    return norm(sub(segment.end,segment.start));
}

double dc_resistance(const FilamentSegment& segment) {
    validate(segment);
    return filament_length(segment)/(segment.conductivity_s_per_m*std::numbers::pi*segment.radius_m*segment.radius_m);
}

double skin_depth_m(double conductivity,double frequency,double relative_permeability) {
    if (!(conductivity>0.0) || !(frequency>0.0) || !(relative_permeability>0.0)
        || !std::isfinite(conductivity) || !std::isfinite(frequency) || !std::isfinite(relative_permeability)) {
        throw std::invalid_argument("invalid skin-depth parameters");
    }
    const double omega=2.0*std::numbers::pi*frequency;
    return std::sqrt(2.0/(omega*mu0*relative_permeability*conductivity));
}

double round_wire_ac_resistance(const FilamentSegment& segment,double frequency,double relative_permeability) {
    validate(segment);
    if (frequency==0.0) return dc_resistance(segment);
    const double delta=skin_depth_m(segment.conductivity_s_per_m,frequency,relative_permeability);
    const double inner=std::max(0.0,segment.radius_m-delta);
    const double effective_area=std::numbers::pi*(segment.radius_m*segment.radius_m-inner*inner);
    return filament_length(segment)/(segment.conductivity_s_per_m*effective_area);
}

double self_partial_inductance(const FilamentSegment& segment) {
    validate(segment);
    const double length=filament_length(segment);
    // High-frequency external partial inductance of a slender round filament.
    return mu0*length/(2.0*std::numbers::pi)*(std::log(2.0*length/segment.radius_m)-1.0);
}

double mutual_partial_inductance(const FilamentSegment& first,const FilamentSegment& second,std::size_t order) {
    validate(first); validate(second);
    const Point3 d1=sub(first.end,first.start), d2=sub(second.end,second.start);
    const double tangent_dot=dot(d1,d2);
    if (std::abs(tangent_dot)<1.0e-30) return 0.0;
    const auto rule=gauss_rule(order);
    double integral=0.0;
    for (std::size_t i=0; i<rule.x.size(); ++i) {
        const double u=0.5*(rule.x[i]+1.0);
        const Point3 p1=add(first.start,scale(d1,u));
        for (std::size_t j=0; j<rule.x.size(); ++j) {
            const double v=0.5*(rule.x[j]+1.0);
            const Point3 p2=add(second.start,scale(d2,v));
            const double distance=norm(sub(p1,p2));
            if (distance<0.25*std::min(first.radius_m,second.radius_m)) {
                throw std::invalid_argument("PEEC mutual integral has overlapping filaments");
            }
            integral+=0.25*rule.w[i]*rule.w[j]/distance;
        }
    }
    return mu0/(4.0*std::numbers::pi)*tangent_dot*integral;
}

void PeecFilamentSystem::add_segment(FilamentSegment segment) {
    validate(segment);
    segments_.push_back(segment);
}

PeecMatrices PeecFilamentSystem::extract(std::size_t order) const {
    if (segments_.empty()) throw std::runtime_error("PEEC extraction requires at least one segment");
    const std::size_t n=segments_.size();
    PeecMatrices out;
    out.size=n;
    out.resistance_ohm.resize(n);
    out.partial_inductance_h.assign(n*n,0.0);
    for (std::size_t i=0; i<n; ++i) {
        out.resistance_ohm[i]=dc_resistance(segments_[i]);
        out.partial_inductance_h[i*n+i]=self_partial_inductance(segments_[i]);
        for (std::size_t j=0; j<i; ++j) {
            const double mutual=mutual_partial_inductance(segments_[i],segments_[j],order);
            out.partial_inductance_h[i*n+j]=out.partial_inductance_h[j*n+i]=mutual;
        }
    }
    return out;
}


PeecCapacitanceMatrices PeecFilamentSystem::extract_capacitance(double relative_permittivity,std::size_t order) const {
    if(segments_.empty())throw std::runtime_error("PEEC capacitance extraction requires at least one segment");
    if(!(relative_permittivity>0.0)||!std::isfinite(relative_permittivity))
        throw std::invalid_argument("PEEC relative permittivity must be finite and positive");
    const std::size_t n=segments_.size();
    PeecCapacitanceMatrices out;out.size=n;out.coefficient_of_potential_v_per_c.assign(n*n,0.0);
    for(std::size_t i=0;i<n;++i){
        out.coefficient_of_potential_v_per_c[i*n+i]=potential_coefficient(segments_[i],segments_[i],relative_permittivity,order,true);
        for(std::size_t j=0;j<i;++j){
            const double value=potential_coefficient(segments_[i],segments_[j],relative_permittivity,order,false);
            out.coefficient_of_potential_v_per_c[i*n+j]=out.coefficient_of_potential_v_per_c[j*n+i]=value;
        }
    }
    out.capacitance_f=inverse_dense_real(out.coefficient_of_potential_v_per_c,n);
    return out;
}

std::vector<std::complex<double>> PeecFilamentSystem::impedance_matrix(double frequency,std::size_t order) const {
    if (!(frequency>=0.0) || !std::isfinite(frequency)) throw std::invalid_argument("invalid PEEC frequency");
    const auto data=extract(order);
    const double omega=2.0*std::numbers::pi*frequency;
    std::vector<std::complex<double>> impedance(data.size*data.size);
    for (std::size_t row=0; row<data.size; ++row) {
        for (std::size_t column=0; column<data.size; ++column) {
            impedance[row*data.size+column]={row==column?data.resistance_ohm[row]:0.0,
                                             omega*data.partial_inductance_h[row*data.size+column]};
        }
    }
    return impedance;
}

std::vector<std::complex<double>> PeecFilamentSystem::impedance_matrix_skin_effect(double frequency,std::size_t order) const {
    if (!(frequency>=0.0) || !std::isfinite(frequency)) throw std::invalid_argument("invalid PEEC frequency");
    const auto data=extract(order);
    const double omega=2.0*std::numbers::pi*frequency;
    std::vector<std::complex<double>> impedance(data.size*data.size);
    for (std::size_t row=0; row<data.size; ++row) {
        for (std::size_t column=0; column<data.size; ++column) {
            const double resistance=row==column?round_wire_ac_resistance(segments_[row],frequency):0.0;
            impedance[row*data.size+column]={resistance,omega*data.partial_inductance_h[row*data.size+column]};
        }
    }
    return impedance;
}

std::string PeecFilamentSystem::to_spice_subcircuit(std::string_view name,std::size_t order) const {
    if (segments_.empty()) throw std::runtime_error("PEEC export requires at least one segment");
    if (name.empty()) throw std::invalid_argument("PEEC subcircuit name must not be empty");
    const auto data=extract(order);
    std::ostringstream out;
    out<<std::setprecision(17)<<".SUBCKT "<<name;
    for (std::size_t i=0;i<data.size;++i) out<<" P"<<i+1U<<" N"<<i+1U;
    out<<'\n';
    for (std::size_t i=0;i<data.size;++i) {
        out<<"R"<<i+1U<<" P"<<i+1U<<" X"<<i+1U<<' '<<data.resistance_ohm[i]<<'\n';
        out<<"L"<<i+1U<<" X"<<i+1U<<" N"<<i+1U<<' '<<data.partial_inductance_h[i*data.size+i]<<'\n';
    }
    for (std::size_t i=0;i<data.size;++i) for (std::size_t j=0;j<i;++j) {
        const double denominator=std::sqrt(data.partial_inductance_h[i*data.size+i]*data.partial_inductance_h[j*data.size+j]);
        if (denominator>0.0) {
            const double coupling=data.partial_inductance_h[i*data.size+j]/denominator;
            if (std::abs(coupling)>1.0e-15) out<<"K"<<j+1U<<'_'<<i+1U<<" L"<<j+1U<<" L"<<i+1U<<' '<<coupling<<'\n';
        }
    }
    out<<".ENDS "<<name<<'\n';
    return out.str();
}

std::vector<std::complex<double>> PeecFilamentSystem::solve_currents(
    double frequency,std::span<const std::complex<double>> voltage,std::size_t order) const {
    if (voltage.size()!=segments_.size()) throw std::invalid_argument("PEEC branch-voltage size mismatch");
    return solve_dense(impedance_matrix(frequency,order),
                       std::vector<std::complex<double>>(voltage.begin(),voltage.end()),segments_.size());
}

} // namespace cfd::rf

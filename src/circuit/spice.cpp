#include "cfd/circuit/spice.hpp"
#include "cfd/circuit/analysis.hpp"
#include "cfd/core/csr_matrix.hpp"
#include "cfd/core/iterative_solvers.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <unordered_set>

namespace cfd::circuit {
namespace {

constexpr double boltzmann = 1.380649e-23;
constexpr double electron_charge = 1.602176634e-19;

std::string upper(std::string text) {
    for (char& c : text) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return text;
}

std::string trim(std::string text) {
    const auto not_space = [](unsigned char c) { return !std::isspace(c); };
    text.erase(text.begin(), std::find_if(text.begin(), text.end(), not_space));
    text.erase(std::find_if(text.rbegin(), text.rend(), not_space).base(), text.end());
    return text;
}

std::vector<std::string> tokens(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string token;
    while (stream >> token) result.push_back(token);
    return result;
}

std::unordered_map<std::string,std::string> model_params(std::string line) {
    for (char& c : line) if (c == '(' || c == ')' || c == ',') c = ' ';
    std::unordered_map<std::string,std::string> result;
    for (const auto& token : tokens(line)) {
        const auto separator = token.find('=');
        if (separator != std::string::npos) result[upper(token.substr(0, separator))] = token.substr(separator + 1U);
    }
    return result;
}

template<class T>
std::vector<T> solve_dense(std::vector<T> matrix, std::vector<T> rhs, std::size_t n) {
    if (matrix.size() != n * n || rhs.size() != n) throw std::invalid_argument("dense circuit matrix size mismatch");
    for (std::size_t column = 0; column < n; ++column) {
        std::size_t pivot = column;
        double best = std::abs(matrix[column * n + column]);
        for (std::size_t row = column + 1U; row < n; ++row) {
            const double candidate = std::abs(matrix[row * n + column]);
            if (candidate > best) { best = candidate; pivot = row; }
        }
        if (best < 1.0e-24) throw std::runtime_error("singular circuit matrix");
        if (pivot != column) {
            for (std::size_t j = column; j < n; ++j) std::swap(matrix[column * n + j], matrix[pivot * n + j]);
            std::swap(rhs[column], rhs[pivot]);
        }
        const T diagonal = matrix[column * n + column];
        for (std::size_t j = column; j < n; ++j) matrix[column * n + j] /= diagonal;
        rhs[column] /= diagonal;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == column) continue;
            const T factor = matrix[row * n + column];
            if (std::abs(factor) == 0.0) continue;
            for (std::size_t j = column; j < n; ++j) matrix[row * n + j] -= factor * matrix[column * n + j];
            rhs[row] -= factor * rhs[column];
        }
    }
    return rhs;
}

class RealMnaMatrix {
public:
    using value_type = double;

    RealMnaMatrix(std::size_t dimension, bool sparse)
        : dimension_(dimension), sparse_(sparse) {
        if (dimension_ == 0U) throw std::invalid_argument("MNA matrix dimension must be positive");
        if (!sparse_) dense_.assign(dimension_ * dimension_, 0.0);
    }

    [[nodiscard]] bool sparse() const noexcept { return sparse_; }
    [[nodiscard]] std::size_t dimension() const noexcept { return dimension_; }

    double& operator[](std::size_t index) {
        if (index >= dimension_ * dimension_) throw std::out_of_range("MNA matrix index out of range");
        if (!sparse_) return dense_[index];
        return sparse_entries_[index];
    }

    [[nodiscard]] std::vector<double> dense_copy() const {
        if (!sparse_) return dense_;
        std::vector<double> result(dimension_ * dimension_, 0.0);
        for (const auto& [index, value] : sparse_entries_) result[index] = value;
        return result;
    }

    [[nodiscard]] cfd::core::CsrMatrix csr() const {
        if (!sparse_) throw std::logic_error("dense MNA matrix has no sparse stamp storage");
        cfd::core::CsrBuilder builder(dimension_, dimension_);
        // Preserve the structural diagonal even for ideal voltage-source branch
        // equations. CsrBuilder intentionally drops exact zeros, while ILU(0)
        // requires a diagonal slot that it can safely floor during factorization.
        const double structural_zero = std::numeric_limits<double>::denorm_min();
        for (std::size_t row = 0; row < dimension_; ++row) builder.add(row, row, structural_zero);
        for (const auto& [index, value] : sparse_entries_) {
            builder.add(index / dimension_, index % dimension_, value);
        }
        return builder.build();
    }

private:
    std::size_t dimension_{};
    bool sparse_{};
    std::vector<double> dense_;
    std::unordered_map<std::size_t, double> sparse_entries_;
};

std::vector<double> solve_real_mna(const RealMnaMatrix& matrix, const std::vector<double>& rhs,
                                   const NewtonConfig& config) {
    const std::size_t n = matrix.dimension();
    if (rhs.size() != n) throw std::invalid_argument("MNA right-hand-side size mismatch");
    if (!matrix.sparse()) return solve_dense(matrix.dense_copy(), rhs, n);
    if (config.sparse_max_iterations == 0U || config.sparse_gmres_restart == 0U
        || !(config.sparse_relative_tolerance > 0.0))
        throw std::invalid_argument("invalid sparse MNA controls");

    try {
        const auto csr = matrix.csr();
        cfd::core::Ilu0Preconditioner ilu(csr);
        std::vector<double> solution(n, 0.0);
        const auto result = cfd::core::restarted_gmres(
            rhs, solution,
            [&](std::span<const double> input, std::span<double> output) { csr.multiply(input, output); },
            [&](std::span<const double> input, std::span<double> output) { ilu(input, output); },
            config.sparse_max_iterations, std::min(config.sparse_gmres_restart, n),
            config.sparse_relative_tolerance);
        if (result.converged) return solution;
        if (!config.sparse_dense_fallback) throw std::runtime_error("sparse MNA GMRES did not converge");
    } catch (const std::exception&) {
        if (!config.sparse_dense_fallback) throw;
    }
    return solve_dense(matrix.dense_copy(), rhs, n);
}

template<class Matrix>
void stamp_admittance(Matrix& matrix, std::size_t dimension,
                      std::size_t positive, std::size_t negative, typename Matrix::value_type admittance) {
    if (positive) matrix[(positive - 1U) * dimension + (positive - 1U)] += admittance;
    if (negative) matrix[(negative - 1U) * dimension + (negative - 1U)] += admittance;
    if (positive && negative) {
        matrix[(positive - 1U) * dimension + (negative - 1U)] -= admittance;
        matrix[(negative - 1U) * dimension + (positive - 1U)] -= admittance;
    }
}

template<class T>
void stamp_current(std::vector<T>& rhs, std::size_t positive, std::size_t negative, T current) {
    if (positive) rhs[positive - 1U] -= current;
    if (negative) rhs[negative - 1U] += current;
}

template<class Matrix>
void stamp_voltage_branch(Matrix& matrix, std::vector<typename Matrix::value_type>& rhs,
                          std::size_t dimension, std::size_t node_count,
                          std::size_t branch, std::size_t positive, std::size_t negative,
                          typename Matrix::value_type voltage,
                          typename Matrix::value_type branch_diagonal = typename Matrix::value_type{}) {
    using T = typename Matrix::value_type;
    const std::size_t equation = (node_count - 1U) + branch;
    if (positive) {
        const std::size_t row = positive - 1U;
        matrix[row * dimension + equation] += T{1};
        matrix[equation * dimension + row] += T{1};
    }
    if (negative) {
        const std::size_t row = negative - 1U;
        matrix[row * dimension + equation] -= T{1};
        matrix[equation * dimension + row] -= T{1};
    }
    matrix[equation * dimension + equation] += branch_diagonal;
    rhs[equation] += voltage;
}

template<class Matrix>
void stamp_vccs(Matrix& matrix, std::size_t dimension,
                std::size_t p, std::size_t n, std::size_t cp, std::size_t cn,
                typename Matrix::value_type gain) {
    using T = typename Matrix::value_type;
    const auto add = [&](std::size_t row, std::size_t column, T value) {
        if (row && column) matrix[(row - 1U) * dimension + (column - 1U)] += value;
    };
    add(p, cp, gain); add(p, cn, -gain);
    add(n, cp, -gain); add(n, cn, gain);
}

template<class Matrix>
void stamp_cccs(Matrix& matrix, std::size_t dimension, std::size_t node_count,
                std::size_t p, std::size_t n, std::size_t control_branch,
                typename Matrix::value_type gain) {
    const std::size_t column = (node_count - 1U) + control_branch;
    if (p) matrix[(p - 1U) * dimension + column] += gain;
    if (n) matrix[(n - 1U) * dimension + column] -= gain;
}

template<class Matrix>
void stamp_vcvs_control(Matrix& matrix, std::size_t dimension, std::size_t node_count,
                        std::size_t branch, std::size_t cp, std::size_t cn,
                        typename Matrix::value_type gain) {
    const std::size_t row = (node_count - 1U) + branch;
    if (cp) matrix[row * dimension + (cp - 1U)] -= gain;
    if (cn) matrix[row * dimension + (cn - 1U)] += gain;
}

template<class Matrix>
void stamp_ccvs_control(Matrix& matrix, std::size_t dimension, std::size_t node_count,
                        std::size_t branch, std::size_t control_branch,
                        typename Matrix::value_type transresistance) {
    const std::size_t row = (node_count - 1U) + branch;
    const std::size_t column = (node_count - 1U) + control_branch;
    matrix[row * dimension + column] -= transresistance;
}

double node_voltage(const std::vector<double>& state, std::size_t node) {
    return node ? state[node - 1U] : 0.0;
}

template<class Matrix>
void stamp_static_device(Matrix& matrix,std::vector<typename Matrix::value_type>* rhs,std::size_t dimension,
                         std::span<const std::size_t> terminals,const StaticDeviceEvaluator& evaluator,
                         std::span<const double> node_voltages) {
    using T = typename Matrix::value_type;
    std::vector<double> local(terminals.size());
    for(std::size_t i=0;i<terminals.size();++i){
        if(terminals[i]>=node_voltages.size())throw std::out_of_range("static-device terminal node out of range");
        local[i]=node_voltages[terminals[i]];
    }
    const auto evaluation=evaluator(local);
    if(evaluation.terminal_current.size()!=terminals.size()
       ||evaluation.jacobian.size()!=terminals.size()*terminals.size())
        throw std::runtime_error("static-device evaluator returned invalid dimensions");
    for(std::size_t r=0;r<terminals.size();++r){
        const auto row_node=terminals[r];if(!row_node)continue;
        double ieq=evaluation.terminal_current[r];
        for(std::size_t c=0;c<terminals.size();++c){
            const double g=evaluation.jacobian[r*terminals.size()+c];
            if(!std::isfinite(g))throw std::runtime_error("static-device Jacobian is non-finite");
            ieq-=g*local[c];
            if(terminals[c])matrix[(row_node-1U)*dimension+(terminals[c]-1U)]+=T{g};
        }
        if(!std::isfinite(ieq))throw std::runtime_error("static-device equivalent current is non-finite");
        if(rhs)(*rhs)[row_node-1U]-=T{ieq};
    }
}


DynamicDeviceEvaluation evaluate_dynamic_device(std::span<const std::size_t> terminals,
                                                const DynamicDeviceEvaluator& evaluator,
                                                std::span<const double> node_voltages) {
    if(!evaluator)throw std::invalid_argument("empty dynamic-device evaluator");
    std::vector<double> local(terminals.size());
    for(std::size_t i=0;i<terminals.size();++i){
        if(terminals[i]>=node_voltages.size())throw std::out_of_range("dynamic-device terminal node out of range");
        local[i]=node_voltages[terminals[i]];
    }
    auto evaluation=evaluator(local);
    const std::size_t n=terminals.size(),matrix_size=n*n;
    if(evaluation.terminal_current.size()!=n||evaluation.terminal_charge.size()!=n
       ||evaluation.current_jacobian.size()!=matrix_size||evaluation.charge_jacobian.size()!=matrix_size)
        throw std::runtime_error("dynamic-device evaluator returned invalid dimensions");
    const auto finite=[](const std::vector<double>& values){return std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);});};
    if(!finite(evaluation.terminal_current)||!finite(evaluation.terminal_charge)
       ||!finite(evaluation.current_jacobian)||!finite(evaluation.charge_jacobian))
        throw std::runtime_error("dynamic-device evaluator returned non-finite values");
    return evaluation;
}

template<class Matrix>
void stamp_dynamic_dc(Matrix& matrix,std::vector<double>& rhs,std::size_t dimension,
                      std::span<const std::size_t> terminals,const DynamicDeviceEvaluator& evaluator,
                      std::span<const double> node_voltages) {
    using T=typename Matrix::value_type;
    const auto evaluation=evaluate_dynamic_device(terminals,evaluator,node_voltages);
    std::vector<double> local(terminals.size());
    for(std::size_t i=0;i<terminals.size();++i)local[i]=node_voltages[terminals[i]];
    for(std::size_t r=0;r<terminals.size();++r){
        const auto row_node=terminals[r];if(!row_node)continue;
        double ieq=evaluation.terminal_current[r];
        for(std::size_t c=0;c<terminals.size();++c){
            const double g=evaluation.current_jacobian[r*terminals.size()+c];
            ieq-=g*local[c];
            if(terminals[c])matrix[(row_node-1U)*dimension+(terminals[c]-1U)]+=T{g};
        }
        rhs[row_node-1U]-=ieq;
    }
}

template<class Matrix>
void stamp_dynamic_ac(Matrix& matrix,std::size_t dimension,
                      std::span<const std::size_t> terminals,const DynamicDeviceEvaluator& evaluator,
                      std::span<const double> node_voltages,double omega) {
    using T=typename Matrix::value_type;
    static_assert(std::is_same_v<T,Complex>,"dynamic AC stamping requires complex matrix");
    const auto evaluation=evaluate_dynamic_device(terminals,evaluator,node_voltages);
    const Complex jw{0.0,omega};
    for(std::size_t r=0;r<terminals.size();++r){
        if(!terminals[r])continue;
        for(std::size_t c=0;c<terminals.size();++c){
            if(!terminals[c])continue;
            const std::size_t k=r*terminals.size()+c;
            matrix[(terminals[r]-1U)*dimension+(terminals[c]-1U)]
                +=Complex{evaluation.current_jacobian[k],0.0}+jw*evaluation.charge_jacobian[k];
        }
    }
}

template<class Matrix>
void stamp_dynamic_transient(Matrix& matrix,std::vector<double>& rhs,std::size_t dimension,
                             std::span<const std::size_t> terminals,const DynamicDeviceEvaluator& evaluator,
                             std::span<const double> current_nodes,std::span<const double> old_nodes,
                             const std::vector<double>* older_nodes,double a0,double a1,double a2) {
    using T=typename Matrix::value_type;
    const auto current=evaluate_dynamic_device(terminals,evaluator,current_nodes);
    const auto old=evaluate_dynamic_device(terminals,evaluator,old_nodes);
    DynamicDeviceEvaluation older;
    if(older_nodes)older=evaluate_dynamic_device(terminals,evaluator,*older_nodes);
    std::vector<double> local(terminals.size());
    for(std::size_t i=0;i<terminals.size();++i)local[i]=current_nodes[terminals[i]];
    for(std::size_t r=0;r<terminals.size();++r){
        const auto row_node=terminals[r];if(!row_node)continue;
        double residual=current.terminal_current[r]+a0*current.terminal_charge[r]+a1*old.terminal_charge[r];
        if(older_nodes)residual+=a2*older.terminal_charge[r];
        double ieq=residual;
        for(std::size_t c=0;c<terminals.size();++c){
            const std::size_t k=r*terminals.size()+c;
            const double jac=current.current_jacobian[k]+a0*current.charge_jacobian[k];
            ieq-=jac*local[c];
            if(terminals[c])matrix[(row_node-1U)*dimension+(terminals[c]-1U)]+=T{jac};
        }
        rhs[row_node-1U]-=ieq;
    }
}

struct BehavioralLinearization {
    double value{};
    std::vector<double> derivative;
};

BehavioralLinearization behavioral_linearize(const BehavioralEvaluator& evaluator,
                                             std::span<const double> node_values,double time) {
    if(!evaluator)throw std::invalid_argument("empty behavioral-source evaluator");
    BehavioralLinearization result;result.value=evaluator(node_values,time);
    if(!std::isfinite(result.value))throw std::runtime_error("behavioral source returned non-finite value");
    result.derivative.assign(node_values.size(),0.0);
    std::vector<double> perturbed(node_values.begin(),node_values.end());
    for(std::size_t node=1U;node<node_values.size();++node){
        const double h=1.0e-7*std::max(1.0,std::abs(node_values[node]));
        perturbed[node]=node_values[node]+h;const double plus=evaluator(perturbed,time);
        perturbed[node]=node_values[node]-h;const double minus=evaluator(perturbed,time);
        perturbed[node]=node_values[node];
        if(!std::isfinite(plus)||!std::isfinite(minus))throw std::runtime_error("behavioral-source Jacobian is non-finite");
        result.derivative[node]=(plus-minus)/(2.0*h);
    }
    return result;
}

template<class Matrix>
void stamp_behavioral_current(Matrix& matrix,std::vector<typename Matrix::value_type>* rhs,std::size_t dimension,
                              std::size_t p,std::size_t n,const BehavioralLinearization& linear,
                              std::span<const double> node_values,double scale=1.0) {
    using T = typename Matrix::value_type;
    double ieq=scale*linear.value;
    for(std::size_t node=1U;node<linear.derivative.size();++node){
        const double g=scale*linear.derivative[node];
        ieq-=g*node_values[node];
        if(p)matrix[(p-1U)*dimension+(node-1U)]+=T{g};
        if(n)matrix[(n-1U)*dimension+(node-1U)]-=T{g};
    }
    if(rhs)stamp_current(*rhs,p,n,T{ieq});
}

template<class Matrix>
void stamp_behavioral_voltage(Matrix& matrix,std::vector<typename Matrix::value_type>* rhs,std::size_t dimension,
                              std::size_t node_count,std::size_t branch,
                              std::size_t p,std::size_t n,const BehavioralLinearization& linear,
                              std::span<const double> node_values,double scale=1.0) {
    using T = typename Matrix::value_type;
    std::vector<T> dummy(dimension,T{});
    auto& target_rhs=rhs?*rhs:dummy;
    stamp_voltage_branch(matrix,target_rhs,dimension,node_count,branch,p,n,T{});
    const std::size_t row=(node_count-1U)+branch;
    double ieq=scale*linear.value;
    for(std::size_t node=1U;node<linear.derivative.size();++node){
        const double g=scale*linear.derivative[node];
        ieq-=g*node_values[node];
        matrix[row*dimension+(node-1U)]-=T{g};
    }
    if(rhs)(*rhs)[row]+=T{ieq};
}

std::vector<double> node_values_from_state(const std::vector<double>& state,std::size_t node_count) {
    std::vector<double> values(node_count,0.0);
    for(std::size_t node=1U;node<node_count;++node)values[node]=state[node-1U];
    return values;
}

struct MosLinearization { double id{}, gm{}, gds{}, ieq{}; };

MosLinearization mos_linearize(const MosLevel1Model& model, double vd, double vg, double vs) {
    const double polarity = model.pmos ? -1.0 : 1.0;
    const double vgs = polarity * (vg - vs);
    const double vds = polarity * (vd - vs);
    const double threshold = std::abs(model.threshold_voltage);
    if (vgs <= threshold || vds <= 0.0) return {};
    const double overdrive = vgs - threshold;
    const double kp = model.transconductance_parameter;
    const double lambda = model.channel_length_modulation;
    double effective_id = 0.0, gm = 0.0, gds = 0.0;
    if (vds < overdrive) {
        const double base = kp * (overdrive * vds - 0.5 * vds * vds);
        effective_id = base * (1.0 + lambda * vds);
        gm = kp * vds * (1.0 + lambda * vds);
        gds = kp * ((overdrive - vds) * (1.0 + lambda * vds)
                    + (overdrive * vds - 0.5 * vds * vds) * lambda);
    } else {
        effective_id = 0.5 * kp * overdrive * overdrive * (1.0 + lambda * vds);
        gm = kp * overdrive * (1.0 + lambda * vds);
        gds = 0.5 * kp * overdrive * overdrive * lambda;
    }
    const double id = polarity * effective_id;
    const double physical_vgs = vg - vs;
    const double physical_vds = vd - vs;
    return {id, gm, gds, id - gm * physical_vgs - gds * physical_vds};
}

MosLinearization jfet_linearize(const JfetModel& model, double vd, double vg, double vs) {
    const double polarity = model.pchannel ? -1.0 : 1.0;
    const double vgs = polarity * (vg - vs);
    const double vds = polarity * (vd - vs);
    const double pinch = model.pchannel ? -std::abs(model.pinch_off_voltage) : model.pinch_off_voltage;
    if (vgs <= pinch || vds <= 0.0) return {};
    const double overdrive = vgs - pinch;
    const double beta = model.beta;
    const double lambda = model.channel_length_modulation;
    double effective_id = 0.0, gm = 0.0, gds = 0.0;
    if (vds < overdrive) {
        const double base = beta * (overdrive * vds - 0.5 * vds * vds);
        effective_id = base * (1.0 + lambda * vds);
        gm = beta * vds * (1.0 + lambda * vds);
        gds = beta * ((overdrive - vds) * (1.0 + lambda * vds)
                      + (overdrive * vds - 0.5 * vds * vds) * lambda);
    } else {
        effective_id = 0.5 * beta * overdrive * overdrive * (1.0 + lambda * vds);
        gm = beta * overdrive * (1.0 + lambda * vds);
        gds = 0.5 * beta * overdrive * overdrive * lambda;
    }
    const double id = polarity * effective_id;
    return {id, gm, gds, id - gm * (vg - vs) - gds * (vd - vs)};
}

template<class Matrix>
void stamp_transistor_linearization(Matrix& matrix, std::size_t dimension,
                                    const MosLinearization& model,
                                    std::size_t drain, std::size_t gate, std::size_t source) {
    using T = typename Matrix::value_type;
    const auto add = [&](std::size_t row, std::size_t column, T value) {
        if (row && column) matrix[(row - 1U) * dimension + (column - 1U)] += value;
    };
    add(drain, drain, T{model.gds});
    add(drain, gate, T{model.gm});
    add(drain, source, T{-model.gds - model.gm});
    add(source, drain, T{-model.gds});
    add(source, gate, T{-model.gm});
    add(source, source, T{model.gds + model.gm});
}

std::array<double,3> bjt_currents(const BjtModel& model, double vc, double vb, double ve) {
    const double polarity = model.pnp ? -1.0 : 1.0;
    const double vt = boltzmann * model.temperature_k / electron_charge;
    const double vbe = polarity * (vb - ve);
    const double vbc = polarity * (vb - vc);
    const double forward = model.saturation_current * (std::exp(std::clamp(vbe / vt, -40.0, 40.0)) - 1.0);
    const double reverse = model.saturation_current * (std::exp(std::clamp(vbc / vt, -40.0, 40.0)) - 1.0);
    const double ic_effective = forward - reverse * (1.0 + 1.0 / model.reverse_beta);
    const double ib_effective = forward / model.forward_beta + reverse / model.reverse_beta;
    const double ie_effective = -ic_effective - ib_effective;
    return {polarity * ic_effective, polarity * ib_effective, polarity * ie_effective};
}

struct BjtLinearization {
    std::array<double,3> ieq{};
    std::array<std::array<double,3>,3> jacobian{};
};

BjtLinearization bjt_linearize(const BjtModel& model, double vc, double vb, double ve) {
    const std::array<double,3> voltage{vc,vb,ve};
    const auto current = bjt_currents(model,vc,vb,ve);
    BjtLinearization result;
    constexpr double h = 1.0e-6;
    for (std::size_t column = 0; column < 3U; ++column) {
        auto plus = voltage, minus = voltage;
        plus[column] += h; minus[column] -= h;
        const auto ip = bjt_currents(model,plus[0],plus[1],plus[2]);
        const auto im = bjt_currents(model,minus[0],minus[1],minus[2]);
        for (std::size_t row=0; row<3U; ++row) result.jacobian[row][column]=(ip[row]-im[row])/(2.0*h);
    }
    for (std::size_t row=0; row<3U; ++row) {
        result.ieq[row]=current[row];
        for (std::size_t column=0; column<3U; ++column) result.ieq[row]-=result.jacobian[row][column]*voltage[column];
    }
    return result;
}

template<class Matrix>
void stamp_bjt_jacobian(Matrix& matrix,std::size_t dimension,const BjtLinearization& model,
                        std::size_t collector,std::size_t base,std::size_t emitter) {
    using T = typename Matrix::value_type;
    const std::array<std::size_t,3> nodes{collector,base,emitter};
    for (std::size_t row=0; row<3U; ++row) {
        if (!nodes[row]) continue;
        for (std::size_t column=0; column<3U; ++column) {
            if (!nodes[column]) continue;
            matrix[(nodes[row]-1U)*dimension+(nodes[column]-1U)] += T{model.jacobian[row][column]};
        }
    }
}

template<class Matrix>
void stamp_bjt_dc(Matrix& matrix,std::vector<double>& rhs,std::size_t dimension,
                  const BjtLinearization& model,std::size_t collector,std::size_t base,std::size_t emitter) {
    stamp_bjt_jacobian(matrix,dimension,model,collector,base,emitter);
    const std::array<std::size_t,3> nodes{collector,base,emitter};
    for (std::size_t terminal=0; terminal<3U; ++terminal)
        if (nodes[terminal]) rhs[nodes[terminal]-1U]-=model.ieq[terminal];
}

double switch_conductance(const SwitchModel& model, double control_voltage) {
    const double upper = model.threshold_voltage + std::abs(model.hysteresis_voltage);
    const double lower = model.threshold_voltage - std::abs(model.hysteresis_voltage);
    const bool on = control_voltage >= upper || (control_voltage > lower && control_voltage >= model.threshold_voltage);
    return 1.0 / (on ? model.on_resistance : model.off_resistance);
}

cfd::rf::NPortPoint interpolate_nport(const std::vector<cfd::rf::NPortPoint>& samples,double frequency) {
    if(samples.empty()) throw std::invalid_argument("N-port requires at least one frequency sample");
    if(samples.size()==1U || frequency<=samples.front().frequency_hz) return samples.front();
    if(frequency>=samples.back().frequency_hz) return samples.back();
    const auto upper_it=std::upper_bound(samples.begin(),samples.end(),frequency,
        [](double value,const cfd::rf::NPortPoint& point){return value<point.frequency_hz;});
    const auto lower_it=upper_it-1;
    if(lower_it->s.size()!=upper_it->s.size() || lower_it->reference_impedance!=upper_it->reference_impedance)
        throw std::runtime_error("N-port interpolation topology/reference mismatch");
    const double fraction=(frequency-lower_it->frequency_hz)/(upper_it->frequency_hz-lower_it->frequency_hz);
    cfd::rf::NPortPoint result;result.frequency_hz=frequency;result.reference_impedance=lower_it->reference_impedance;result.s=cfd::rf::ComplexMatrix(lower_it->s.size());
    for(std::size_t row=0;row<result.s.size();++row)for(std::size_t column=0;column<result.s.size();++column)
        result.s(row,column)=lower_it->s(row,column)*(1.0-fraction)+upper_it->s(row,column)*fraction;
    return result;
}

template<class Device>
void stamp_nport(std::vector<Complex>& matrix,std::size_t dimension,const Device& device,double frequency) {
    const auto point=interpolate_nport(device.samples,frequency);
    if(point.s.size()!=device.ports.size()) throw std::runtime_error("N-port device port count mismatch");
    const auto z=cfd::rf::s_to_z(point.s,point.reference_impedance);
    const auto y=cfd::rf::inverse(z);
    const auto add=[&](std::size_t row,std::size_t column,Complex value){if(row&&column)matrix[(row-1U)*dimension+(column-1U)]+=value;};
    for(std::size_t i=0;i<device.ports.size();++i){
        const auto [pi,ni]=device.ports[i];
        for(std::size_t j=0;j<device.ports.size();++j){
            const auto [pj,nj]=device.ports[j];const Complex value=y(i,j);
            add(pi,pj,value);add(pi,nj,-value);add(ni,pj,-value);add(ni,nj,value);
        }
    }
}

} // namespace

double PulseWaveform::operator()(double time) const {
    if (!std::isfinite(time) || !std::isfinite(initial) || !std::isfinite(pulsed)
        || !(delay >= 0.0) || !(rise >= 0.0) || !(fall >= 0.0) || !(width >= 0.0) || !(period >= 0.0)) {
        throw std::invalid_argument("invalid PULSE waveform");
    }
    if (time < delay) return initial;
    double local = time - delay;
    if (period > 0.0) local = std::fmod(local, period);
    if (rise > 0.0 && local < rise) return initial + (pulsed - initial) * (local / rise);
    local -= rise;
    if (local < width) return pulsed;
    local -= width;
    if (fall > 0.0 && local < fall) return pulsed + (initial - pulsed) * (local / fall);
    return initial;
}

double SineWaveform::operator()(double time) const {
    if (!std::isfinite(time) || !std::isfinite(offset) || !std::isfinite(amplitude)
        || !(frequency_hz >= 0.0) || !(delay >= 0.0) || !(damping_per_s >= 0.0) || !std::isfinite(phase_deg)) {
        throw std::invalid_argument("invalid SIN waveform");
    }
    if (time < delay) return offset;
    const double local = time - delay;
    const double phase = phase_deg * std::numbers::pi / 180.0;
    return offset + amplitude * std::exp(-damping_per_s * local)
        * std::sin(2.0 * std::numbers::pi * frequency_hz * local + phase);
}

double PwlWaveform::operator()(double time) const {
    if (!std::isfinite(time) || points.empty()) throw std::invalid_argument("invalid PWL waveform");
    if (time <= points.front().first) return points.front().second;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const auto [t0,v0] = points[i-1U];
        const auto [t1,v1] = points[i];
        if (!(t1 > t0) || !std::isfinite(t0) || !std::isfinite(t1) || !std::isfinite(v0) || !std::isfinite(v1)) {
            throw std::invalid_argument("PWL times must be finite and strictly increasing");
        }
        if (time <= t1) return v0 + (v1-v0) * ((time-t0)/(t1-t0));
    }
    return points.back().second;
}

Circuit::Circuit() {
    node_names_.push_back("0");
    node_lookup_["0"] = 0U;
    node_lookup_["GND"] = 0U;
}

std::size_t Circuit::node(std::string_view name) {
    std::string key = upper(std::string(name));
    if (key == "GND") key = "0";
    if (const auto found = node_lookup_.find(key); found != node_lookup_.end()) return found->second;
    const std::size_t id = node_names_.size();
    node_names_.push_back(std::string(name));
    node_lookup_[key] = id;
    return id;
}

std::size_t Circuit::find_node(std::string_view name) const {
    std::string key = upper(std::string(name));
    if (key == "GND") key = "0";
    const auto found = node_lookup_.find(key);
    if (found == node_lookup_.end()) throw std::out_of_range("unknown circuit node");
    return found->second;
}

void Circuit::add_resistor(std::string name,std::size_t p,std::size_t n,double value) {
    if (!(value > 0.0) || !std::isfinite(value)) throw std::invalid_argument("invalid resistance");
    resistors_.push_back({std::move(name),p,n,value});
}
void Circuit::add_capacitor(std::string name,std::size_t p,std::size_t n,double value) {
    if (!(value >= 0.0) || !std::isfinite(value)) throw std::invalid_argument("invalid capacitance");
    capacitors_.push_back({std::move(name),p,n,value});
}
void Circuit::add_inductor(std::string name,std::size_t p,std::size_t n,double value) {
    if (!(value > 0.0) || !std::isfinite(value)) throw std::invalid_argument("invalid inductance");
    inductors_.push_back({std::move(name),p,n,value});
}
void Circuit::add_mutual_inductance(std::string name,std::string first,std::string second,double coupling) {
    if (!std::isfinite(coupling) || std::abs(coupling) > 1.0) throw std::invalid_argument("mutual-inductor coupling must be in [-1,1]");
    mutual_inductances_.push_back({std::move(name),upper(std::move(first)),upper(std::move(second)),coupling});
}
void Circuit::add_voltage_source(std::string name,std::size_t p,std::size_t n,double dc,Complex ac_value,SourceWaveform transient_waveform) {
    if (!std::isfinite(dc) || !std::isfinite(ac_value.real()) || !std::isfinite(ac_value.imag())) throw std::invalid_argument("invalid voltage source");
    voltage_sources_.push_back({std::move(name),p,n,dc,ac_value,std::move(transient_waveform)});
}
void Circuit::add_current_source(std::string name,std::size_t p,std::size_t n,double dc,Complex ac_value,SourceWaveform transient_waveform) {
    if (!std::isfinite(dc) || !std::isfinite(ac_value.real()) || !std::isfinite(ac_value.imag())) throw std::invalid_argument("invalid current source");
    current_sources_.push_back({std::move(name),p,n,dc,ac_value,std::move(transient_waveform)});
}
void Circuit::add_vccs(std::string name,std::size_t p,std::size_t n,std::size_t cp,std::size_t cn,double gain) {
    if (!std::isfinite(gain)) throw std::invalid_argument("invalid VCCS gain");
    vccs_.push_back({std::move(name),p,n,cp,cn,gain});
}
void Circuit::add_vcvs(std::string name,std::size_t p,std::size_t n,std::size_t cp,std::size_t cn,double gain) {
    if (!std::isfinite(gain)) throw std::invalid_argument("invalid VCVS gain");
    vcvs_.push_back({std::move(name),p,n,cp,cn,gain});
}
void Circuit::add_cccs(std::string name,std::size_t p,std::size_t n,std::string control,double gain) {
    if (!std::isfinite(gain)) throw std::invalid_argument("invalid CCCS gain");
    cccs_.push_back({std::move(name),p,n,upper(std::move(control)),gain});
}
void Circuit::add_ccvs(std::string name,std::size_t p,std::size_t n,std::string control,double gain) {
    if (!std::isfinite(gain)) throw std::invalid_argument("invalid CCVS transresistance");
    ccvs_.push_back({std::move(name),p,n,upper(std::move(control)),gain});
}
void Circuit::add_ideal_transformer(std::string name,std::size_t pp,std::size_t pn,
                                    std::size_t sp,std::size_t sn,double turns_ratio) {
    if(!(turns_ratio>0.0)||!std::isfinite(turns_ratio))throw std::invalid_argument("invalid ideal-transformer turns ratio");
    const std::string voltage_branch=name+"$V";
    add_vcvs(voltage_branch,pp,pn,sp,sn,turns_ratio);
    add_cccs(name+"$I",sp,sn,voltage_branch,-turns_ratio);
}
void Circuit::add_diode_model(std::string name,DiodeModel model) {
    if (!(model.saturation_current>0.0) || !(model.emission_coefficient>0.0) || !(model.temperature_k>0.0)) throw std::invalid_argument("invalid diode model");
    diode_models_[upper(std::move(name))]=model;
}
void Circuit::add_diode(std::string name,std::size_t a,std::size_t k,std::string model) { diodes_.push_back({std::move(name),a,k,upper(std::move(model))}); }
void Circuit::add_mos_model(std::string name,MosLevel1Model model) {
    if (!(model.transconductance_parameter>0.0) || model.channel_length_modulation<0.0 || !std::isfinite(model.threshold_voltage)) throw std::invalid_argument("invalid MOS level-1 model");
    mos_models_[upper(std::move(name))]=model;
}
void Circuit::add_mosfet(std::string name,std::size_t d,std::size_t g,std::size_t s,std::size_t b,std::string model) { mosfets_.push_back({std::move(name),d,g,s,b,upper(std::move(model))}); }
void Circuit::add_bjt_model(std::string name,BjtModel model) {
    if (!(model.saturation_current>0.0) || !(model.forward_beta>0.0) || !(model.reverse_beta>0.0) || !(model.temperature_k>0.0)) throw std::invalid_argument("invalid BJT model");
    bjt_models_[upper(std::move(name))]=model;
}
void Circuit::add_bjt(std::string name,std::size_t c,std::size_t b,std::size_t e,std::string model) { bjts_.push_back({std::move(name),c,b,e,upper(std::move(model))}); }
void Circuit::add_jfet_model(std::string name,JfetModel model) {
    if (!(model.beta>0.0) || model.channel_length_modulation<0.0 || !std::isfinite(model.pinch_off_voltage)) throw std::invalid_argument("invalid JFET model");
    jfet_models_[upper(std::move(name))]=model;
}
void Circuit::add_jfet(std::string name,std::size_t d,std::size_t g,std::size_t s,std::string model) { jfets_.push_back({std::move(name),d,g,s,upper(std::move(model))}); }
void Circuit::add_switch_model(std::string name,SwitchModel model) {
    if (!(model.on_resistance>0.0) || !(model.off_resistance>model.on_resistance) || !std::isfinite(model.threshold_voltage) || !std::isfinite(model.hysteresis_voltage)) throw std::invalid_argument("invalid voltage-controlled switch model");
    switch_models_[upper(std::move(name))]=model;
}
void Circuit::add_switch(std::string name,std::size_t p,std::size_t n,std::size_t cp,std::size_t cn,std::string model) { switches_.push_back({std::move(name),p,n,cp,cn,upper(std::move(model))}); }
void Circuit::add_nport(std::string name,std::vector<std::pair<std::size_t,std::size_t>> ports,std::vector<cfd::rf::NPortPoint> samples) {
    if(ports.empty()||samples.empty()) throw std::invalid_argument("N-port device requires ports and samples");
    std::sort(samples.begin(),samples.end(),[](const auto& a,const auto& b){return a.frequency_hz<b.frequency_hz;});
    for(const auto& point:samples){if(!(point.frequency_hz>0.0)||point.s.size()!=ports.size()||point.reference_impedance.size()!=ports.size())throw std::invalid_argument("invalid N-port sample");}
    nports_.push_back({std::move(name),std::move(ports),std::move(samples)});
}
void Circuit::add_tem_transmission_line(std::string name,
                                         std::pair<std::size_t,std::size_t> first_port,
                                         std::pair<std::size_t,std::size_t> second_port,
                                         const cfd::rf::TemTransmissionLine& line,
                                         std::span<const double> frequencies,
                                         double reference_impedance) {
    if (frequencies.empty() || !(reference_impedance>0.0) || !std::isfinite(reference_impedance))
        throw std::invalid_argument("transmission-line samples/reference impedance are invalid");
    std::vector<cfd::rf::NPortPoint> samples;
    samples.reserve(frequencies.size());
    for (double frequency : frequencies) {
        const auto two=line.s_parameters(frequency,reference_impedance);
        cfd::rf::NPortPoint point;
        point.frequency_hz=frequency;
        point.reference_impedance={reference_impedance,reference_impedance};
        point.s=cfd::rf::ComplexMatrix(2U);
        point.s(0,0)=two.a11;point.s(0,1)=two.a12;
        point.s(1,0)=two.a21;point.s(1,1)=two.a22;
        samples.push_back(std::move(point));
    }
    add_nport(std::move(name),{first_port,second_port},std::move(samples));
}

void Circuit::add_microstrip_line(std::string name,
                                   std::pair<std::size_t,std::size_t> first_port,
                                   std::pair<std::size_t,std::size_t> second_port,
                                   double width,double height,double er,double length,
                                   std::span<const double> frequencies,double reference_impedance) {
    if (frequencies.empty() || !(length>=0.0) || !std::isfinite(length))
        throw std::invalid_argument("microstrip circuit element requires samples and non-negative length");
    const auto quasi=cfd::rf::microstrip_quasi_static(width,height,er,frequencies.front());
    cfd::rf::TemTransmissionLine line{quasi.characteristic_impedance,quasi.phase_velocity,0.0,length};
    add_tem_transmission_line(std::move(name),first_port,second_port,line,frequencies,reference_impedance);
}

void Circuit::add_stripline_line(std::string name,
                                 std::pair<std::size_t,std::size_t> first_port,
                                 std::pair<std::size_t,std::size_t> second_port,
                                 double width,double spacing,double er,double length,
                                 std::span<const double> frequencies,double reference_impedance) {
    if(frequencies.empty()||!(length>=0.0)||!std::isfinite(length))
        throw std::invalid_argument("stripline circuit element requires samples and non-negative length");
    const auto quasi=cfd::rf::stripline_quasi_static(width,spacing,er,frequencies.front());
    cfd::rf::TemTransmissionLine line{quasi.characteristic_impedance,quasi.phase_velocity,0.0,length};
    add_tem_transmission_line(std::move(name),first_port,second_port,line,frequencies,reference_impedance);
}

void Circuit::add_coplanar_waveguide_line(std::string name,
                                          std::pair<std::size_t,std::size_t> first_port,
                                          std::pair<std::size_t,std::size_t> second_port,
                                          double width,double gap,double er,double length,
                                          std::span<const double> frequencies,double reference_impedance) {
    if(frequencies.empty()||!(length>=0.0)||!std::isfinite(length))
        throw std::invalid_argument("coplanar-waveguide circuit element requires samples and non-negative length");
    const auto quasi=cfd::rf::coplanar_waveguide_quasi_static(width,gap,er,frequencies.front());
    cfd::rf::TemTransmissionLine line{quasi.characteristic_impedance,quasi.phase_velocity,0.0,length};
    add_tem_transmission_line(std::move(name),first_port,second_port,line,frequencies,reference_impedance);
}

void Circuit::add_rectangular_waveguide_te10_line(std::string name,
                                                   std::pair<std::size_t,std::size_t> first_port,
                                                   std::pair<std::size_t,std::size_t> second_port,
                                                   double broad,double narrow,double er,double mur,double length,
                                                   std::span<const double> frequencies,double reference_impedance) {
    if(frequencies.empty()||!(length>=0.0)||!std::isfinite(length)||!(reference_impedance>0.0))
        throw std::invalid_argument("waveguide circuit element requires samples and non-negative length");
    std::vector<cfd::rf::NPortPoint> samples;
    samples.reserve(frequencies.size());
    for(double frequency:frequencies){
        const auto mode=cfd::rf::rectangular_waveguide_te10(broad,narrow,er,mur,frequency);
        cfd::rf::TransmissionLine line{mode.wave_impedance,{0.0,mode.propagation_constant_rad_per_m},length};
        const auto two=line.s_parameters(reference_impedance);
        cfd::rf::NPortPoint point;point.frequency_hz=frequency;point.s=cfd::rf::ComplexMatrix(2U);
        point.reference_impedance={reference_impedance,reference_impedance};
        point.s(0,0)=two.a11;point.s(0,1)=two.a12;point.s(1,0)=two.a21;point.s(1,1)=two.a22;
        samples.push_back(std::move(point));
    }
    add_nport(std::move(name),{first_port,second_port},std::move(samples));
}
void Circuit::add_static_device(std::string name,std::vector<std::size_t> terminals,StaticDeviceEvaluator evaluator) {
    if(terminals.empty()||!evaluator)throw std::invalid_argument("static device requires terminals and evaluator");
    for(const auto node_id:terminals)if(node_id>=node_count())throw std::out_of_range("static-device terminal node out of range");
    static_devices_.push_back({std::move(name),std::move(terminals),std::move(evaluator)});
}

void Circuit::add_dynamic_device(std::string name,std::vector<std::size_t> terminals,DynamicDeviceEvaluator evaluator) {
    if(terminals.empty()||!evaluator)throw std::invalid_argument("dynamic device requires terminals and evaluator");
    for(const auto node_id:terminals)if(node_id>=node_count())throw std::out_of_range("dynamic-device terminal node out of range");
    dynamic_devices_.push_back({std::move(name),std::move(terminals),std::move(evaluator)});
}

void Circuit::add_electrothermal_device(std::string name,std::vector<std::size_t> electrical_terminals,
                                        std::size_t thermal_node,ElectroThermalDeviceConfig config,
                                        ElectroThermalDeviceEvaluator evaluator) {
    if(name.empty()||electrical_terminals.empty()||!evaluator||thermal_node==0U||thermal_node>=node_count()
       ||!(config.ambient_temperature_k>0.0)||!std::isfinite(config.ambient_temperature_k)
       ||!(config.thermal_resistance_k_per_w>0.0)||!std::isfinite(config.thermal_resistance_k_per_w)
       ||config.thermal_capacitance_j_per_k<0.0||!std::isfinite(config.thermal_capacitance_j_per_k))
        throw std::invalid_argument("invalid electrothermal device configuration");
    for(const auto node_id:electrical_terminals)if(node_id>=node_count())throw std::out_of_range("electrothermal terminal node out of range");
    std::vector<std::size_t> dae_terminals=electrical_terminals;dae_terminals.push_back(thermal_node);
    const std::size_t electrical_count=electrical_terminals.size();
    add_dynamic_device(std::move(name),std::move(dae_terminals),
        [electrical_count,config,evaluator=std::move(evaluator)](std::span<const double> local){
            if(local.size()!=electrical_count+1U)throw std::runtime_error("electrothermal local state size mismatch");
            const double theta=local[electrical_count];
            const double temperature=config.ambient_temperature_k+theta;
            if(!(temperature>0.0)||!std::isfinite(temperature))throw std::runtime_error("electrothermal device reached non-physical temperature");
            const auto electrical=evaluator(local.first(electrical_count),temperature);
            if(electrical.electrical_current.size()!=electrical_count
               ||electrical.electrical_jacobian.size()!=electrical_count*electrical_count
               ||electrical.electrical_temperature_derivative.size()!=electrical_count
               ||electrical.power_voltage_derivative.size()!=electrical_count
               ||!std::isfinite(electrical.dissipated_power_w)
               ||!std::isfinite(electrical.power_temperature_derivative_w_per_k))
                throw std::runtime_error("electrothermal evaluator returned invalid dimensions/values");
            DynamicDeviceEvaluation result;
            const std::size_t total=electrical_count+1U;
            result.terminal_current.assign(total,0.0);result.current_jacobian.assign(total*total,0.0);
            result.terminal_charge.assign(total,0.0);result.charge_jacobian.assign(total*total,0.0);
            for(std::size_t r=0;r<electrical_count;++r){
                if(!std::isfinite(electrical.electrical_current[r])||!std::isfinite(electrical.electrical_temperature_derivative[r]))
                    throw std::runtime_error("electrothermal electrical current is non-finite");
                result.terminal_current[r]=electrical.electrical_current[r];
                for(std::size_t c=0;c<electrical_count;++c){
                    const double g=electrical.electrical_jacobian[r*electrical_count+c];
                    if(!std::isfinite(g))throw std::runtime_error("electrothermal electrical Jacobian is non-finite");
                    result.current_jacobian[r*total+c]=g;
                }
                result.current_jacobian[r*total+electrical_count]=electrical.electrical_temperature_derivative[r];
            }
            const std::size_t thermal=electrical_count;
            result.terminal_current[thermal]=theta/config.thermal_resistance_k_per_w-electrical.dissipated_power_w;
            for(std::size_t c=0;c<electrical_count;++c){
                const double derivative=electrical.power_voltage_derivative[c];
                if(!std::isfinite(derivative))throw std::runtime_error("electrothermal power Jacobian is non-finite");
                result.current_jacobian[thermal*total+c]=-derivative;
            }
            result.current_jacobian[thermal*total+thermal]
                =1.0/config.thermal_resistance_k_per_w-electrical.power_temperature_derivative_w_per_k;
            result.terminal_charge[thermal]=config.thermal_capacitance_j_per_k*theta;
            result.charge_jacobian[thermal*total+thermal]=config.thermal_capacitance_j_per_k;
            return result;
        });
}

void Circuit::add_behavioral_current_source(std::string name,std::size_t positive,std::size_t negative,BehavioralEvaluator evaluator) {
    if(name.empty()||!evaluator||positive>=node_count()||negative>=node_count())
        throw std::invalid_argument("invalid behavioral current source");
    behavioral_current_sources_.push_back({std::move(name),positive,negative,std::move(evaluator)});
}

void Circuit::add_behavioral_voltage_source(std::string name,std::size_t positive,std::size_t negative,BehavioralEvaluator evaluator) {
    if(name.empty()||!evaluator||positive>=node_count()||negative>=node_count())
        throw std::invalid_argument("invalid behavioral voltage source");
    behavioral_voltage_sources_.push_back({std::move(name),positive,negative,std::move(evaluator)});
}

void Circuit::set_voltage_source_dc(std::string_view name,double voltage) {
    if (!std::isfinite(voltage)) throw std::invalid_argument("invalid voltage-source DC value");
    const std::string target=upper(std::string(name));
    for (auto& source:voltage_sources_) if (upper(source.name)==target) { source.dc=voltage; return; }
    throw std::out_of_range("unknown voltage source");
}

double Circuit::voltage_source_dc(std::string_view name) const {
    const std::string target=upper(std::string(name));
    for(const auto& source:voltage_sources_)if(upper(source.name)==target)return source.dc;
    throw std::out_of_range("unknown voltage source");
}
void Circuit::set_resistance(std::string_view name,double resistance_value) {
    if(!(resistance_value>0.0)||!std::isfinite(resistance_value))throw std::invalid_argument("invalid resistance");
    const std::string target=upper(std::string(name));for(auto& resistor:resistors_)if(upper(resistor.name)==target){resistor.value=resistance_value;return;}throw std::out_of_range("unknown resistor");
}
double Circuit::resistance(std::string_view name) const {
    const std::string target=upper(std::string(name));for(const auto& resistor:resistors_)if(upper(resistor.name)==target)return resistor.value;throw std::out_of_range("unknown resistor");
}

void Circuit::set_device_temperature(double temperature) {
    if (!(temperature > 0.0) || !std::isfinite(temperature))
        throw std::invalid_argument("device temperature must be finite and positive");
    for (auto& [name, model] : diode_models_) { (void)name; model.temperature_k = temperature; }
    for (auto& [name, model] : bjt_models_) { (void)name; model.temperature_k = temperature; }
}

void Circuit::set_initial_voltage(std::size_t node_id,double voltage) {
    if(node_id>=node_count()||!std::isfinite(voltage))throw std::invalid_argument("invalid circuit initial condition");
    if(node_id==0U){if(voltage!=0.0)throw std::invalid_argument("ground initial condition must be zero");return;}
    initial_voltage_[node_id]=voltage;
}

void Circuit::set_nodeset_voltage(std::size_t node_id,double voltage) {
    if(node_id>=node_count()||!std::isfinite(voltage))throw std::invalid_argument("invalid circuit nodeset");
    if(node_id==0U){if(voltage!=0.0)throw std::invalid_argument("ground nodeset must be zero");return;}
    nodeset_voltage_[node_id]=voltage;
}


std::size_t Circuit::branch_count() const noexcept { return voltage_sources_.size()+behavioral_voltage_sources_.size()+vcvs_.size()+ccvs_.size()+inductors_.size(); }
std::size_t Circuit::inductor_branch(std::size_t index) const noexcept { return voltage_sources_.size()+behavioral_voltage_sources_.size()+vcvs_.size()+ccvs_.size()+index; }

std::size_t Circuit::branch_index(std::string_view name) const {
    const std::string target=upper(std::string(name));
    std::size_t branch=0U;
    for (const auto& source:voltage_sources_) { if (upper(source.name)==target) return branch; ++branch; }
    for (const auto& source:behavioral_voltage_sources_) { if (upper(source.name)==target) return branch; ++branch; }
    for (const auto& source:vcvs_) { if (upper(source.name)==target) return branch; ++branch; }
    for (const auto& source:ccvs_) { if (upper(source.name)==target) return branch; ++branch; }
    for (const auto& inductor:inductors_) { if (upper(inductor.name)==target) return branch; ++branch; }
    throw std::out_of_range("unknown branch-current control element");
}

OperatingPoint Circuit::dc_operating_point_scaled(const NewtonConfig& cfg,double source_scale,const OperatingPoint* initial) const {
    if (cfg.max_iterations==0U || !(cfg.voltage_tolerance>0.0) || !(cfg.residual_tolerance>0.0) || cfg.gmin<0.0 || !std::isfinite(source_scale)) throw std::invalid_argument("invalid circuit Newton controls");
    const std::size_t node_unknowns=node_count()-1U;
    const std::size_t dimension=node_unknowns+branch_count();
    std::vector<double> state(dimension,0.0);
    if (initial && initial->node_voltage.size()==node_count() && initial->branch_current.size()==branch_count()) {
        for (std::size_t node_id=1; node_id<node_count(); ++node_id) state[node_id-1U]=initial->node_voltage[node_id];
        std::copy(initial->branch_current.begin(),initial->branch_current.end(),state.begin()+static_cast<std::ptrdiff_t>(node_unknowns));
    }
    else {
        for(const auto& [node_id,voltage]:nodeset_voltage_)
            if(node_id>0U&&node_id<node_count())state[node_id-1U]=voltage;
    }
    std::vector<double> previous=state;

    for (std::size_t iteration=0; iteration<cfg.max_iterations; ++iteration) {
        const bool use_sparse=cfg.sparse_mna_threshold>0U&&dimension>=cfg.sparse_mna_threshold;
        RealMnaMatrix matrix(dimension,use_sparse);std::vector<double> rhs(dimension,0.0);
        for (const auto& resistor:resistors_) stamp_admittance(matrix,dimension,resistor.p,resistor.n,1.0/resistor.value);
        for (std::size_t node_id=1; node_id<node_count(); ++node_id) matrix[(node_id-1U)*dimension+(node_id-1U)]+=cfg.gmin;
        for (const auto& source:current_sources_) stamp_current(rhs,source.p,source.n,source_scale*source.dc);
        for (const auto& source:vccs_) stamp_vccs(matrix,dimension,source.p,source.n,source.cp,source.cn,source.gain);
        const auto behavioral_nodes=node_values_from_state(previous,node_count());
        for(const auto& source:behavioral_current_sources_){
            const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,0.0);
            stamp_behavioral_current(matrix,&rhs,dimension,source.p,source.n,linear,behavioral_nodes,source_scale);
        }

        std::size_t branch=0U;
        for (const auto& source:voltage_sources_) {
            stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,source_scale*source.dc);
            ++branch;
        }
        for(const auto& source:behavioral_voltage_sources_){
            const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,0.0);
            stamp_behavioral_voltage(matrix,&rhs,dimension,node_count(),branch,source.p,source.n,linear,behavioral_nodes,source_scale);
            ++branch;
        }
        for (const auto& source:vcvs_) {
            stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,0.0);
            stamp_vcvs_control(matrix,dimension,node_count(),branch,source.cp,source.cn,source.gain);
            ++branch;
        }
        for (const auto& source:ccvs_) {
            stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,0.0);
            stamp_ccvs_control(matrix,dimension,node_count(),branch,branch_index(source.control),source.gain);
            ++branch;
        }
        for (const auto& inductor:inductors_) {
            stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,inductor.p,inductor.n,0.0);
            ++branch;
        }
        for (const auto& source:cccs_) stamp_cccs(matrix,dimension,node_count(),source.p,source.n,branch_index(source.control),source.gain);

        for (const auto& diode:diodes_) {
            const auto found=diode_models_.find(diode.model); if(found==diode_models_.end())throw std::runtime_error("unknown diode model: "+diode.model);
            const auto& model=found->second;
            const double vt=model.emission_coefficient*boltzmann*model.temperature_k/electron_charge;
            const double vd=node_voltage(previous,diode.a)-node_voltage(previous,diode.k);
            const double exponential=std::exp(std::clamp(vd/vt,-40.0,40.0));
            const double current=model.saturation_current*(exponential-1.0);
            const double conductance=model.saturation_current*exponential/vt;
            stamp_admittance(matrix,dimension,diode.a,diode.k,conductance);
            stamp_current(rhs,diode.a,diode.k,current-conductance*vd);
        }
        for (const auto& device:mosfets_) {
            const auto found=mos_models_.find(device.model); if(found==mos_models_.end())throw std::runtime_error("unknown MOS model: "+device.model);
            const auto linear=mos_linearize(found->second,node_voltage(previous,device.d),node_voltage(previous,device.g),node_voltage(previous,device.s));
            stamp_transistor_linearization(matrix,dimension,linear,device.d,device.g,device.s);
            stamp_current(rhs,device.d,device.s,linear.ieq);
        }
        for (const auto& device:jfets_) {
            const auto found=jfet_models_.find(device.model); if(found==jfet_models_.end())throw std::runtime_error("unknown JFET model: "+device.model);
            const auto linear=jfet_linearize(found->second,node_voltage(previous,device.d),node_voltage(previous,device.g),node_voltage(previous,device.s));
            stamp_transistor_linearization(matrix,dimension,linear,device.d,device.g,device.s);
            stamp_current(rhs,device.d,device.s,linear.ieq);
        }
        for (const auto& device:bjts_) {
            const auto found=bjt_models_.find(device.model); if(found==bjt_models_.end())throw std::runtime_error("unknown BJT model: "+device.model);
            const auto linear=bjt_linearize(found->second,node_voltage(previous,device.c),node_voltage(previous,device.b),node_voltage(previous,device.e));
            stamp_bjt_dc(matrix,rhs,dimension,linear,device.c,device.b,device.e);
        }
        for (const auto& device:switches_) {
            const auto found=switch_models_.find(device.model); if(found==switch_models_.end())throw std::runtime_error("unknown switch model: "+device.model);
            const double control=node_voltage(previous,device.cp)-node_voltage(previous,device.cn);
            stamp_admittance(matrix,dimension,device.p,device.n,switch_conductance(found->second,control));
        }
        if(!static_devices_.empty()){
            std::vector<double> node_values(node_count(),0.0);
            for(std::size_t node_id=1;node_id<node_count();++node_id)node_values[node_id]=previous[node_id-1U];
            for(const auto& device:static_devices_)stamp_static_device(matrix,&rhs,dimension,device.terminals,device.evaluator,node_values);
        }
        if(!dynamic_devices_.empty()){
            const auto node_values=node_values_from_state(previous,node_count());
            for(const auto& device:dynamic_devices_)stamp_dynamic_dc(matrix,rhs,dimension,device.terminals,device.evaluator,node_values);
        }

        auto candidate=solve_real_mna(matrix,rhs,cfg);
        // Limit PN-junction Newton voltage jumps before the global damping step.
        for(const auto& diode:diodes_){
            const auto& model=diode_models_.at(diode.model);
            const double vt=model.emission_coefficient*boltzmann*model.temperature_k/electron_charge;
            const double old_v=node_voltage(previous,diode.a)-node_voltage(previous,diode.k);
            const double proposed=node_voltage(candidate,diode.a)-node_voltage(candidate,diode.k);
            const double limited=limit_pn_junction_voltage(proposed,old_v,vt,model.saturation_current);
            const double correction=limited-proposed;
            if(diode.a&&diode.k){candidate[diode.a-1U]+=0.5*correction;candidate[diode.k-1U]-=0.5*correction;}
            else if(diode.a)candidate[diode.a-1U]+=correction;
            else if(diode.k)candidate[diode.k-1U]-=correction;
        }
        // Nonlinear device state is controlled by node voltages. Branch-current
        // unknowns (voltage sources/inductors) can legitimately be very large in
        // low-impedance circuits and must not throttle the Newton voltage step.
        double maximum_voltage_delta=0.0;
        for (std::size_t i=0;i<node_unknowns;++i)
            maximum_voltage_delta=std::max(maximum_voltage_delta,std::abs(candidate[i]-previous[i]));
        const bool nonlinear = !diodes_.empty() || !mosfets_.empty() || !jfets_.empty() ||
            !bjts_.empty() || !switches_.empty() || !static_devices_.empty() || !dynamic_devices_.empty() ||
            !behavioral_current_sources_.empty() || !behavioral_voltage_sources_.empty();
        const double damping = nonlinear && maximum_voltage_delta > 0.5
            ? 0.5 / maximum_voltage_delta : 1.0;
        for (std::size_t i=0;i<dimension;++i) state[i]=previous[i]+damping*(candidate[i]-previous[i]);
        previous=state;
        if (maximum_voltage_delta*damping<cfg.voltage_tolerance) {
            OperatingPoint result;
            result.node_voltage.assign(node_count(),0.0);
            for (std::size_t node_id=1;node_id<node_count();++node_id) result.node_voltage[node_id]=state[node_id-1U];
            result.branch_current.assign(state.begin()+static_cast<std::ptrdiff_t>(node_unknowns),state.end());
            result.iterations=iteration+1U; result.converged=true;
            return result;
        }
    }
    OperatingPoint result;
    result.node_voltage.assign(node_count(),0.0);
    for (std::size_t node_id=1;node_id<node_count();++node_id) result.node_voltage[node_id]=previous[node_id-1U];
    result.branch_current.assign(previous.begin()+static_cast<std::ptrdiff_t>(node_unknowns),previous.end());
    result.iterations=cfg.max_iterations; result.converged=false;
    return result;
}

OperatingPoint Circuit::dc_operating_point(const NewtonConfig& cfg) const { return dc_operating_point_scaled(cfg,1.0,nullptr); }

OperatingPoint Circuit::dc_operating_point_homotopy(const NewtonConfig& config,std::size_t source_steps,std::size_t gmin_steps) const {
    if (source_steps==0U || gmin_steps==0U) throw std::invalid_argument("homotopy step counts must be positive");
    NewtonConfig staged=config;
    staged.gmin=std::max(config.gmin,1.0e-3);
    OperatingPoint current;
    for (std::size_t step=0;step<source_steps;++step) {
        const double scale=static_cast<double>(step+1U)/static_cast<double>(source_steps);
        current=dc_operating_point_scaled(staged,scale,current.node_voltage.empty()?nullptr:&current);
        if(!current.converged) return current;
    }
    const double target=std::max(config.gmin,1.0e-15);
    for (std::size_t step=0;step<gmin_steps;++step) {
        const double fraction=static_cast<double>(step+1U)/static_cast<double>(gmin_steps);
        staged.gmin=std::exp(std::log(1.0e-3)*(1.0-fraction)+std::log(target)*fraction);
        current=dc_operating_point_scaled(staged,1.0,&current);
        if(!current.converged) return current;
    }
    return dc_operating_point_scaled(config,1.0,&current);
}

AcSolution Circuit::ac(double frequency,const OperatingPoint* supplied,const NewtonConfig& cfg) const {
    if (!(frequency>0.0) || !std::isfinite(frequency)) throw std::invalid_argument("AC frequency must be positive");
    OperatingPoint owned;
    const OperatingPoint* operating=supplied;
    if (!operating) {
        const bool nonlinear = !diodes_.empty() || !mosfets_.empty() || !bjts_.empty() || !jfets_.empty() || !switches_.empty() || !static_devices_.empty() || !dynamic_devices_.empty() || !behavioral_current_sources_.empty() || !behavioral_voltage_sources_.empty();
        if (nonlinear) {
            owned=dc_operating_point_homotopy(cfg);
            if(!owned.converged) throw std::runtime_error("AC operating point did not converge");
        } else {
            owned.node_voltage.assign(node_count(),0.0);
            owned.branch_current.assign(branch_count(),0.0);
            owned.converged=true;
        }
        operating=&owned;
    }
    const std::size_t node_unknowns=node_count()-1U,dimension=node_unknowns+branch_count();
    std::vector<Complex> matrix(dimension*dimension),rhs(dimension);
    const Complex jw{0.0,2.0*std::numbers::pi*frequency};
    for (const auto& resistor:resistors_) stamp_admittance(matrix,dimension,resistor.p,resistor.n,Complex{1.0/resistor.value,0.0});
    for (const auto& capacitor:capacitors_) stamp_admittance(matrix,dimension,capacitor.p,capacitor.n,jw*capacitor.value);
    for (const auto& source:current_sources_) stamp_current(rhs,source.p,source.n,source.ac);
    for (const auto& source:vccs_) stamp_vccs(matrix,dimension,source.p,source.n,source.cp,source.cn,Complex{source.gain,0.0});
    for(const auto& source:behavioral_current_sources_){
        const auto linear=behavioral_linearize(source.evaluator,operating->node_voltage,0.0);
        stamp_behavioral_current(matrix,static_cast<std::vector<Complex>*>(nullptr),dimension,source.p,source.n,linear,operating->node_voltage);
    }
    std::size_t branch=0U;
    for (const auto& source:voltage_sources_) { stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,source.ac); ++branch; }
    for(const auto& source:behavioral_voltage_sources_){
        const auto linear=behavioral_linearize(source.evaluator,operating->node_voltage,0.0);
        stamp_behavioral_voltage(matrix,static_cast<std::vector<Complex>*>(nullptr),dimension,node_count(),branch,source.p,source.n,linear,operating->node_voltage);
        ++branch;
    }
    for (const auto& source:vcvs_) { stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,Complex{}); stamp_vcvs_control(matrix,dimension,node_count(),branch,source.cp,source.cn,Complex{source.gain,0.0}); ++branch; }
    for (const auto& source:ccvs_) { stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,Complex{}); stamp_ccvs_control(matrix,dimension,node_count(),branch,branch_index(source.control),Complex{source.gain,0.0}); ++branch; }
    for (const auto& inductor:inductors_) { stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,inductor.p,inductor.n,Complex{},-jw*inductor.value); ++branch; }
    for (const auto& coupling:mutual_inductances_) {
        std::size_t first=inductors_.size(),second=inductors_.size();
        for(std::size_t i=0;i<inductors_.size();++i){const std::string name=upper(inductors_[i].name);if(name==coupling.first)first=i;if(name==coupling.second)second=i;}
        if(first==inductors_.size()||second==inductors_.size())throw std::runtime_error("mutual inductance references unknown inductor");
        const double mutual=coupling.coupling*std::sqrt(inductors_[first].value*inductors_[second].value);
        const std::size_t r1=node_unknowns+inductor_branch(first),r2=node_unknowns+inductor_branch(second);
        matrix[r1*dimension+r2]+=-jw*mutual; matrix[r2*dimension+r1]+=-jw*mutual;
    }
    for (const auto& source:cccs_) stamp_cccs(matrix,dimension,node_count(),source.p,source.n,branch_index(source.control),Complex{source.gain,0.0});
    for (const auto& diode:diodes_) {
        const auto& model=diode_models_.at(diode.model);const double vt=model.emission_coefficient*boltzmann*model.temperature_k/electron_charge;
        const double vd=operating->node_voltage[diode.a]-operating->node_voltage[diode.k];
        const double conductance=model.saturation_current*std::exp(std::clamp(vd/vt,-40.0,40.0))/vt;
        stamp_admittance(matrix,dimension,diode.a,diode.k,Complex{conductance,0.0});
    }
    for (const auto& device:mosfets_) {const auto linear=mos_linearize(mos_models_.at(device.model),operating->node_voltage[device.d],operating->node_voltage[device.g],operating->node_voltage[device.s]);stamp_transistor_linearization(matrix,dimension,linear,device.d,device.g,device.s);}
    for (const auto& device:jfets_) {const auto linear=jfet_linearize(jfet_models_.at(device.model),operating->node_voltage[device.d],operating->node_voltage[device.g],operating->node_voltage[device.s]);stamp_transistor_linearization(matrix,dimension,linear,device.d,device.g,device.s);}
    for (const auto& device:bjts_) {const auto linear=bjt_linearize(bjt_models_.at(device.model),operating->node_voltage[device.c],operating->node_voltage[device.b],operating->node_voltage[device.e]);stamp_bjt_jacobian(matrix,dimension,linear,device.c,device.b,device.e);}
    for (const auto& device:switches_) {const double control=operating->node_voltage[device.cp]-operating->node_voltage[device.cn];stamp_admittance(matrix,dimension,device.p,device.n,Complex{switch_conductance(switch_models_.at(device.model),control),0.0});}
    for(const auto& device:static_devices_)stamp_static_device(matrix,static_cast<std::vector<Complex>*>(nullptr),dimension,device.terminals,device.evaluator,operating->node_voltage);
    for(const auto& device:dynamic_devices_)stamp_dynamic_ac(matrix,dimension,device.terminals,device.evaluator,operating->node_voltage,2.0*std::numbers::pi*frequency);
    for (const auto& device:nports_) stamp_nport(matrix,dimension,device,frequency);
    const auto solution=solve_dense(matrix,rhs,dimension);
    AcSolution result;result.frequency_hz=frequency;result.node_voltage.assign(node_count(),{});
    for(std::size_t node_id=1;node_id<node_count();++node_id)result.node_voltage[node_id]=solution[node_id-1U];
    result.branch_current.assign(solution.begin()+static_cast<std::ptrdiff_t>(node_unknowns),solution.end());return result;
}

std::vector<AcSolution> Circuit::ac_log_sweep(double start,double stop,std::size_t points,const NewtonConfig& config) const {
    if (!(start>0.0) || !(stop>=start) || points==0U) throw std::invalid_argument("invalid AC sweep range");
    const auto operating=dc_operating_point_homotopy(config);if(!operating.converged)throw std::runtime_error("AC sweep operating point did not converge");
    std::vector<AcSolution> result;result.reserve(points);
    if(points==1U){result.push_back(ac(start,&operating,config));return result;}
    const double ratio=std::pow(stop/start,1.0/static_cast<double>(points-1U));double frequency=start;
    for(std::size_t i=0;i<points;++i){result.push_back(ac(i+1U==points?stop:frequency,&operating,config));frequency*=ratio;}
    return result;
}

std::vector<DcSweepPoint> Circuit::dc_sweep_voltage_source(std::string_view source,double start,double stop,std::size_t points,const NewtonConfig& config) {
    if(points==0U||!std::isfinite(start)||!std::isfinite(stop))throw std::invalid_argument("invalid DC sweep range");
    const std::string target=upper(std::string(source));VoltageSource* selected=nullptr;
    for(auto& candidate:voltage_sources_)if(upper(candidate.name)==target){selected=&candidate;break;}
    if(!selected) throw std::out_of_range("unknown voltage source for DC sweep");
    const double original=selected->dc;
    std::vector<DcSweepPoint> result;result.reserve(points);
    try{
        for(std::size_t i=0;i<points;++i){const double t=points==1U?0.0:static_cast<double>(i)/static_cast<double>(points-1U);selected->dc=start+(stop-start)*t;result.push_back({selected->dc,dc_operating_point_homotopy(config)});}
    }catch(...){selected->dc=original;throw;}
    selected->dc=original;return result;
}

std::vector<TemperatureSweepPoint> Circuit::dc_temperature_sweep(double start,double stop,std::size_t points,
                                                                  const NewtonConfig& config) const {
    if (points==0U || !(start>0.0) || !(stop>0.0) || !std::isfinite(start) || !std::isfinite(stop))
        throw std::invalid_argument("invalid temperature sweep range");
    std::vector<TemperatureSweepPoint> result;
    result.reserve(points);
    for (std::size_t i=0;i<points;++i) {
        const double fraction=points==1U?0.0:static_cast<double>(i)/static_cast<double>(points-1U);
        const double temperature=start+(stop-start)*fraction;
        Circuit copy=*this;
        copy.set_device_temperature(temperature);
        auto operating=copy.dc_operating_point_homotopy(config);
        result.push_back({temperature,std::move(operating)});
    }
    return result;
}

std::vector<TransientPoint> Circuit::transient(double dt,std::size_t steps,const NewtonConfig& cfg) const {
    return transient(dt,steps,TransientMethod::backward_euler,cfg);
}

std::vector<TransientPoint> Circuit::transient(double dt,std::size_t steps,TransientMethod method,
                                                const NewtonConfig& cfg) const {
    if (!(dt>0.0) || !std::isfinite(dt)) throw std::invalid_argument("transient time step must be positive");
    if(method==TransientMethod::trapezoidal&&!dynamic_devices_.empty())
        throw std::invalid_argument("charge-based dynamic devices currently support backward Euler or BDF2 transient integration");
    const std::size_t node_unknowns=node_count()-1U,dimension=node_unknowns+branch_count();
    std::vector<double> state(dimension,0.0);
    for(const auto& [node_id,voltage]:initial_voltage_)if(node_id>0U&&node_id<node_count())state[node_id-1U]=voltage;
    std::vector<double> older=state,previous=state;
    std::vector<double> capacitor_current(capacitors_.size(),0.0);
    std::vector<double> inductor_voltage(inductors_.size(),0.0);
    std::vector<TransientPoint> result;result.reserve(steps+1U);
    const auto save=[&](double time){
        TransientPoint point;point.time=time;point.node_voltage.assign(node_count(),0.0);
        for(std::size_t node_id=1;node_id<node_count();++node_id)point.node_voltage[node_id]=state[node_id-1U];
        point.branch_current.assign(state.begin()+static_cast<std::ptrdiff_t>(node_unknowns),state.end());
        result.push_back(std::move(point));
    };
    save(0.0);
    for(std::size_t step=0;step<steps;++step){
        const double next_time = static_cast<double>(step + 1U) * dt;
        const auto old=state;
        previous=state;
        const TransientMethod active=(method==TransientMethod::bdf2&&step==0U)
            ?TransientMethod::backward_euler:method;
        bool converged=false;
        for(std::size_t iteration=0;iteration<cfg.max_iterations;++iteration){
            const bool use_sparse=cfg.sparse_mna_threshold>0U&&dimension>=cfg.sparse_mna_threshold;
            RealMnaMatrix matrix(dimension,use_sparse);std::vector<double> rhs(dimension,0.0);
            for(const auto&r:resistors_)stamp_admittance(matrix,dimension,r.p,r.n,1.0/r.value);
            for(std::size_t node_id=1;node_id<node_count();++node_id)
                matrix[(node_id-1U)*dimension+(node_id-1U)]+=cfg.gmin;
            for(const auto& source:current_sources_) {
                const double value = source.transient ? source.transient(next_time) : source.dc;
                stamp_current(rhs,source.p,source.n,value);
            }
            for(const auto&s:vccs_)stamp_vccs(matrix,dimension,s.p,s.n,s.cp,s.cn,s.gain);
            const auto behavioral_nodes=node_values_from_state(previous,node_count());
            for(const auto& source:behavioral_current_sources_){
                const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,next_time);
                stamp_behavioral_current(matrix,&rhs,dimension,source.p,source.n,linear,behavioral_nodes);
            }

            for(std::size_t k=0;k<capacitors_.size();++k){
                const auto& c=capacitors_[k];
                const double v1=node_voltage(old,c.p)-node_voltage(old,c.n);
                double g=0.0,history=0.0;
                if(active==TransientMethod::backward_euler){
                    g=c.value/dt;history=-g*v1;
                }else if(active==TransientMethod::trapezoidal){
                    g=2.0*c.value/dt;history=-g*v1-capacitor_current[k];
                }else{
                    const double v2=node_voltage(older,c.p)-node_voltage(older,c.n);
                    g=1.5*c.value/dt;history=c.value/(2.0*dt)*(-4.0*v1+v2);
                }
                stamp_admittance(matrix,dimension,c.p,c.n,g);
                stamp_current(rhs,c.p,c.n,history);
            }

            std::size_t branch=0U;
            for(const auto& source:voltage_sources_){
                const double value = source.transient ? source.transient(next_time) : source.dc;
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,value);
                ++branch;
            }
            for(const auto& source:behavioral_voltage_sources_){
                const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,next_time);
                stamp_behavioral_voltage(matrix,&rhs,dimension,node_count(),branch,source.p,source.n,linear,behavioral_nodes);
                ++branch;
            }
            for(const auto&e:vcvs_){
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,e.p,e.n,0.0);
                stamp_vcvs_control(matrix,dimension,node_count(),branch,e.cp,e.cn,e.gain);++branch;
            }
            for(const auto&h:ccvs_){
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,h.p,h.n,0.0);
                stamp_ccvs_control(matrix,dimension,node_count(),branch,branch_index(h.control),h.gain);++branch;
            }
            for(std::size_t k=0;k<inductors_.size();++k){
                const auto& l=inductors_[k];const std::size_t ib=inductor_branch(k);
                const double i1=old[node_unknowns+ib];double coefficient=0.0,history=0.0;
                if(active==TransientMethod::backward_euler){
                    coefficient=l.value/dt;history=-coefficient*i1;
                }else if(active==TransientMethod::trapezoidal){
                    coefficient=2.0*l.value/dt;history=-coefficient*i1-inductor_voltage[k];
                }else{
                    const double i2=older[node_unknowns+ib];coefficient=1.5*l.value/dt;
                    history=l.value/(2.0*dt)*(-4.0*i1+i2);
                }
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),ib,l.p,l.n,history,-coefficient);
            }
            for(const auto&coupling:mutual_inductances_){
                std::size_t first=inductors_.size(),second=inductors_.size();
                for(std::size_t i=0;i<inductors_.size();++i){
                    const auto name=upper(inductors_[i].name);if(name==coupling.first)first=i;if(name==coupling.second)second=i;
                }
                if(first==inductors_.size()||second==inductors_.size())
                    throw std::runtime_error("mutual inductance references unknown inductor");
                const double m=coupling.coupling*std::sqrt(inductors_[first].value*inductors_[second].value);
                const std::size_t b1=inductor_branch(first),b2=inductor_branch(second),r1=node_unknowns+b1,r2=node_unknowns+b2;
                double coefficient=0.0,h1=0.0,h2=0.0;
                if(active==TransientMethod::backward_euler){
                    coefficient=m/dt;h1=-coefficient*old[node_unknowns+b2];h2=-coefficient*old[node_unknowns+b1];
                }else if(active==TransientMethod::trapezoidal){
                    coefficient=2.0*m/dt;h1=-coefficient*old[node_unknowns+b2];h2=-coefficient*old[node_unknowns+b1];
                }else{
                    coefficient=1.5*m/dt;
                    h1=m/(2.0*dt)*(-4.0*old[node_unknowns+b2]+older[node_unknowns+b2]);
                    h2=m/(2.0*dt)*(-4.0*old[node_unknowns+b1]+older[node_unknowns+b1]);
                }
                matrix[r1*dimension+r2]-=coefficient;matrix[r2*dimension+r1]-=coefficient;
                rhs[r1]+=h1;rhs[r2]+=h2;
            }
            for(const auto&f:cccs_)stamp_cccs(matrix,dimension,node_count(),f.p,f.n,branch_index(f.control),f.gain);
            for(const auto&d:diodes_){
                const auto&m=diode_models_.at(d.model);const double vt=m.emission_coefficient*boltzmann*m.temperature_k/electron_charge;
                const double raw=node_voltage(previous,d.a)-node_voltage(previous,d.k);
                const double vd=limit_pn_junction_voltage(raw,node_voltage(old,d.a)-node_voltage(old,d.k),vt,m.saturation_current);
                const double ex=std::exp(std::clamp(vd/vt,-40.0,40.0)),current=m.saturation_current*(ex-1.0),g=m.saturation_current*ex/vt;
                stamp_admittance(matrix,dimension,d.a,d.k,g);stamp_current(rhs,d.a,d.k,current-g*vd);
            }
            for(const auto&m:mosfets_){const auto linear=mos_linearize(mos_models_.at(m.model),node_voltage(previous,m.d),node_voltage(previous,m.g),node_voltage(previous,m.s));stamp_transistor_linearization(matrix,dimension,linear,m.d,m.g,m.s);stamp_current(rhs,m.d,m.s,linear.ieq);}
            for(const auto&j:jfets_){const auto linear=jfet_linearize(jfet_models_.at(j.model),node_voltage(previous,j.d),node_voltage(previous,j.g),node_voltage(previous,j.s));stamp_transistor_linearization(matrix,dimension,linear,j.d,j.g,j.s);stamp_current(rhs,j.d,j.s,linear.ieq);}
            for(const auto&q:bjts_){const auto linear=bjt_linearize(bjt_models_.at(q.model),node_voltage(previous,q.c),node_voltage(previous,q.b),node_voltage(previous,q.e));stamp_bjt_dc(matrix,rhs,dimension,linear,q.c,q.b,q.e);}
            for(const auto&s:switches_){const double control=node_voltage(previous,s.cp)-node_voltage(previous,s.cn);stamp_admittance(matrix,dimension,s.p,s.n,switch_conductance(switch_models_.at(s.model),control));}
            if(!static_devices_.empty()){std::vector<double> node_values(node_count(),0.0);for(std::size_t node_id=1;node_id<node_count();++node_id)node_values[node_id]=previous[node_id-1U];for(const auto& device:static_devices_)stamp_static_device(matrix,&rhs,dimension,device.terminals,device.evaluator,node_values);}
            if(!dynamic_devices_.empty()){
                const auto current_nodes=node_values_from_state(previous,node_count());
                const auto old_nodes=node_values_from_state(old,node_count());
                const auto older_nodes=node_values_from_state(older,node_count());
                double a0=1.0/dt,a1=-1.0/dt,a2=0.0;
                if(active==TransientMethod::bdf2){a0=1.5/dt;a1=-2.0/dt;a2=0.5/dt;}
                for(const auto& device:dynamic_devices_)
                    stamp_dynamic_transient(matrix,rhs,dimension,device.terminals,device.evaluator,current_nodes,old_nodes,
                                            active==TransientMethod::bdf2?&older_nodes:nullptr,a0,a1,a2);
            }
            auto candidate=solve_real_mna(matrix,rhs,cfg);double delta=0.0;
            for(std::size_t i=0;i<dimension;++i)delta=std::max(delta,std::abs(candidate[i]-previous[i]));
            const double damping=delta>0.5?0.5/delta:1.0;
            for(std::size_t i=0;i<dimension;++i)state[i]=previous[i]+damping*(candidate[i]-previous[i]);
            previous=state;if(delta*damping<cfg.voltage_tolerance){converged=true;break;}
        }
        if(!converged) throw std::runtime_error("transient circuit Newton solve did not converge");
        if(active==TransientMethod::trapezoidal){
            for(std::size_t k=0;k<capacitors_.size();++k){const auto& c=capacitors_[k];const double v0=node_voltage(old,c.p)-node_voltage(old,c.n),v1=node_voltage(state,c.p)-node_voltage(state,c.n);capacitor_current[k]=2.0*c.value/dt*(v1-v0)-capacitor_current[k];}
            for(std::size_t k=0;k<inductors_.size();++k){const auto& l=inductors_[k];inductor_voltage[k]=node_voltage(state,l.p)-node_voltage(state,l.n);}
        }
        older=old;
        save(next_time);
    }
    return result;
}

AdaptiveTransientResult Circuit::transient_adaptive(double stop_time,
                                                      const AdaptiveTransientConfig& controls,
                                                      const NewtonConfig& cfg) const {
    if (!(stop_time>0.0) || !std::isfinite(stop_time)
        || !(controls.initial_step>0.0) || !(controls.minimum_step>0.0)
        || !(controls.maximum_step>=controls.minimum_step)
        || !(controls.relative_tolerance>0.0) || !(controls.absolute_tolerance>0.0)
        || !(controls.safety_factor>0.0) || !(controls.minimum_scale>0.0)
        || !(controls.maximum_scale>=controls.minimum_scale) || controls.max_step_attempts==0U)
        throw std::invalid_argument("invalid adaptive transient controls");
    if(controls.method==TransientMethod::trapezoidal)
        throw std::invalid_argument("adaptive transient currently supports backward Euler or BDF2");

    const std::size_t node_unknowns=node_count()-1U,dimension=node_unknowns+branch_count();
    std::vector<double> state(dimension,0.0);
    for(const auto& [node_id,voltage]:initial_voltage_)if(node_id>0U&&node_id<node_count())state[node_id-1U]=voltage;
    std::vector<double> older_state=state;

    const auto advance = [&](const std::vector<double>& old,
                             const std::vector<double>* older,
                             double previous_dt,double next_time,double dt,
                             TransientMethod active) {
        double a0=1.0/dt,a1=-1.0/dt,a2=0.0;
        if(active==TransientMethod::bdf2){
            if(older==nullptr || !(previous_dt>0.0) || !std::isfinite(previous_dt))
                throw std::invalid_argument("BDF2 adaptive step requires valid history");
            const double sum=dt+previous_dt;
            a0=(2.0*dt+previous_dt)/(dt*sum);
            a1=-sum/(dt*previous_dt);
            a2=dt/(previous_dt*sum);
        }
        std::vector<double> current=old,previous=old;
        bool converged=false;
        for(std::size_t iteration=0;iteration<cfg.max_iterations;++iteration){
            const bool use_sparse=cfg.sparse_mna_threshold>0U&&dimension>=cfg.sparse_mna_threshold;
            RealMnaMatrix matrix(dimension,use_sparse);std::vector<double> rhs(dimension,0.0);
            for(const auto&r:resistors_)stamp_admittance(matrix,dimension,r.p,r.n,1.0/r.value);
            for(std::size_t node_id=1;node_id<node_count();++node_id)
                matrix[(node_id-1U)*dimension+(node_id-1U)]+=cfg.gmin;
            for(const auto& source:current_sources_){
                const double value=source.transient?source.transient(next_time):source.dc;
                stamp_current(rhs,source.p,source.n,value);
            }
            for(const auto&s:vccs_)stamp_vccs(matrix,dimension,s.p,s.n,s.cp,s.cn,s.gain);
            const auto behavioral_nodes=node_values_from_state(previous,node_count());
            for(const auto& source:behavioral_current_sources_){
                const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,next_time);
                stamp_behavioral_current(matrix,&rhs,dimension,source.p,source.n,linear,behavioral_nodes);
            }
            for(const auto& c:capacitors_){
                const double v1=node_voltage(old,c.p)-node_voltage(old,c.n);
                const double v2=older?(node_voltage(*older,c.p)-node_voltage(*older,c.n)):0.0;
                const double g=c.value*a0;
                const double history=c.value*(a1*v1+a2*v2);
                stamp_admittance(matrix,dimension,c.p,c.n,g);
                stamp_current(rhs,c.p,c.n,history);
            }

            std::size_t branch=0U;
            for(const auto& source:voltage_sources_){
                const double value=source.transient?source.transient(next_time):source.dc;
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,source.p,source.n,value);++branch;
            }
            for(const auto& source:behavioral_voltage_sources_){
                const auto linear=behavioral_linearize(source.evaluator,behavioral_nodes,next_time);
                stamp_behavioral_voltage(matrix,&rhs,dimension,node_count(),branch,source.p,source.n,linear,behavioral_nodes);
                ++branch;
            }
            for(const auto&e:vcvs_){
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,e.p,e.n,0.0);
                stamp_vcvs_control(matrix,dimension,node_count(),branch,e.cp,e.cn,e.gain);++branch;
            }
            for(const auto&h:ccvs_){
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),branch,h.p,h.n,0.0);
                stamp_ccvs_control(matrix,dimension,node_count(),branch,branch_index(h.control),h.gain);++branch;
            }
            for(std::size_t k=0;k<inductors_.size();++k){
                const auto& l=inductors_[k];const std::size_t ib=inductor_branch(k);
                const double i1=old[node_unknowns+ib];
                const double i2=older?(*older)[node_unknowns+ib]:0.0;
                const double coefficient=l.value*a0;
                const double history=l.value*(a1*i1+a2*i2);
                stamp_voltage_branch(matrix,rhs,dimension,node_count(),ib,l.p,l.n,history,-coefficient);
            }
            for(const auto&coupling:mutual_inductances_){
                std::size_t first=inductors_.size(),second=inductors_.size();
                for(std::size_t i=0;i<inductors_.size();++i){
                    const auto name=upper(inductors_[i].name);
                    if(name==coupling.first)first=i;
                    if(name==coupling.second)second=i;
                }
                if(first==inductors_.size()||second==inductors_.size())
                    throw std::runtime_error("mutual inductance references unknown inductor");
                const double mutual=coupling.coupling*std::sqrt(inductors_[first].value*inductors_[second].value);
                const std::size_t b1=inductor_branch(first),b2=inductor_branch(second);
                const std::size_t r1=node_unknowns+b1,r2=node_unknowns+b2;
                const double coefficient=mutual*a0;
                const double old_second=old[node_unknowns+b2],old_first=old[node_unknowns+b1];
                const double older_second=older?(*older)[node_unknowns+b2]:0.0;
                const double older_first=older?(*older)[node_unknowns+b1]:0.0;
                matrix[r1*dimension+r2]-=coefficient;matrix[r2*dimension+r1]-=coefficient;
                rhs[r1]+=mutual*(a1*old_second+a2*older_second);
                rhs[r2]+=mutual*(a1*old_first+a2*older_first);
            }
            for(const auto&f:cccs_)stamp_cccs(matrix,dimension,node_count(),f.p,f.n,branch_index(f.control),f.gain);
            for(const auto&d:diodes_){
                const auto&m=diode_models_.at(d.model);
                const double vt=m.emission_coefficient*boltzmann*m.temperature_k/electron_charge;
                const double raw=node_voltage(previous,d.a)-node_voltage(previous,d.k);
                const double vd=limit_pn_junction_voltage(raw,node_voltage(old,d.a)-node_voltage(old,d.k),vt,m.saturation_current);
                const double ex=std::exp(std::clamp(vd/vt,-40.0,40.0));
                const double diode_current=m.saturation_current*(ex-1.0),g=m.saturation_current*ex/vt;
                stamp_admittance(matrix,dimension,d.a,d.k,g);stamp_current(rhs,d.a,d.k,diode_current-g*vd);
            }
            for(const auto&m:mosfets_){const auto linear=mos_linearize(mos_models_.at(m.model),node_voltage(previous,m.d),node_voltage(previous,m.g),node_voltage(previous,m.s));stamp_transistor_linearization(matrix,dimension,linear,m.d,m.g,m.s);stamp_current(rhs,m.d,m.s,linear.ieq);}
            for(const auto&j:jfets_){const auto linear=jfet_linearize(jfet_models_.at(j.model),node_voltage(previous,j.d),node_voltage(previous,j.g),node_voltage(previous,j.s));stamp_transistor_linearization(matrix,dimension,linear,j.d,j.g,j.s);stamp_current(rhs,j.d,j.s,linear.ieq);}
            for(const auto&q:bjts_){const auto linear=bjt_linearize(bjt_models_.at(q.model),node_voltage(previous,q.c),node_voltage(previous,q.b),node_voltage(previous,q.e));stamp_bjt_dc(matrix,rhs,dimension,linear,q.c,q.b,q.e);}
            for(const auto&s:switches_){const double control=node_voltage(previous,s.cp)-node_voltage(previous,s.cn);stamp_admittance(matrix,dimension,s.p,s.n,switch_conductance(switch_models_.at(s.model),control));}
            if(!static_devices_.empty()){
                std::vector<double> node_values(node_count(),0.0);
                for(std::size_t node_id=1;node_id<node_count();++node_id)node_values[node_id]=previous[node_id-1U];
                for(const auto& device:static_devices_)stamp_static_device(matrix,&rhs,dimension,device.terminals,device.evaluator,node_values);
            }
            if(!dynamic_devices_.empty()){
                const auto current_nodes=node_values_from_state(previous,node_count());
                const auto old_nodes=node_values_from_state(old,node_count());
                std::vector<double> older_nodes;
                if(older)older_nodes=node_values_from_state(*older,node_count());
                for(const auto& device:dynamic_devices_)
                    stamp_dynamic_transient(matrix,rhs,dimension,device.terminals,device.evaluator,current_nodes,old_nodes,
                                            older?&older_nodes:nullptr,a0,a1,a2);
            }
            auto candidate=solve_real_mna(matrix,rhs,cfg);
            double maximum_voltage_delta=0.0;
            for(std::size_t i=0;i<node_unknowns;++i)
                maximum_voltage_delta=std::max(maximum_voltage_delta,std::abs(candidate[i]-previous[i]));
            const bool nonlinear=!diodes_.empty()||!mosfets_.empty()||!jfets_.empty()||!bjts_.empty()||!switches_.empty()||!static_devices_.empty()||!dynamic_devices_.empty()
                ||!behavioral_current_sources_.empty()||!behavioral_voltage_sources_.empty();
            const double damping=nonlinear&&maximum_voltage_delta>0.5?0.5/maximum_voltage_delta:1.0;
            for(std::size_t i=0;i<dimension;++i)current[i]=previous[i]+damping*(candidate[i]-previous[i]);
            previous=current;
            if(maximum_voltage_delta*damping<cfg.voltage_tolerance){converged=true;break;}
        }
        if(!converged)throw std::runtime_error("adaptive transient Newton solve did not converge");
        return current;
    };

    AdaptiveTransientResult result;
    result.minimum_accepted_step=std::numeric_limits<double>::infinity();
    const auto save=[&](double sample_time,const std::vector<double>& accepted){
        TransientPoint point;point.time=sample_time;point.node_voltage.assign(node_count(),0.0);
        for(std::size_t node_id=1;node_id<node_count();++node_id)point.node_voltage[node_id]=accepted[node_id-1U];
        point.branch_current.assign(accepted.begin()+static_cast<std::ptrdiff_t>(node_unknowns),accepted.end());
        result.points.push_back(std::move(point));
    };
    save(0.0,state);
    double time=0.0,previous_dt=0.0;
    bool has_history=false;
    double dt=std::clamp(controls.initial_step,controls.minimum_step,controls.maximum_step);
    std::size_t attempts=0U;
    while(time<stop_time){
        if(++attempts>controls.max_step_attempts)throw std::runtime_error("adaptive transient step-attempt budget exhausted");
        const double remaining=stop_time-time;
        dt=std::min(dt,remaining);
        const bool use_bdf2=controls.method==TransientMethod::bdf2&&has_history;
        const TransientMethod active=use_bdf2?TransientMethod::bdf2:TransientMethod::backward_euler;
        const auto* history=use_bdf2?&older_state:nullptr;
        const double history_dt=use_bdf2?previous_dt:0.0;
        const auto full=advance(state,history,history_dt,time+dt,dt,active);
        const auto half=advance(state,history,history_dt,time+0.5*dt,0.5*dt,active);
        const auto refined=advance(half,use_bdf2?&state:nullptr,use_bdf2?0.5*dt:0.0,
                                   time+dt,0.5*dt,active);
        const double richardson=use_bdf2?3.0:1.0;
        double normalized_error=0.0;
        const std::size_t error_unknowns=node_unknowns?node_unknowns:dimension;
        for(std::size_t i=0;i<error_unknowns;++i){
            const double scale=controls.absolute_tolerance+controls.relative_tolerance*std::max(std::abs(full[i]),std::abs(refined[i]));
            normalized_error=std::max(normalized_error,std::abs(refined[i]-full[i])/(richardson*scale));
        }
        if(!std::isfinite(normalized_error))throw std::runtime_error("adaptive transient produced non-finite LTE estimate");
        if(normalized_error<=1.0){
            older_state=state;state=refined;previous_dt=dt;has_history=true;time+=dt;++result.accepted_steps;
            result.minimum_accepted_step=std::min(result.minimum_accepted_step,dt);
            result.maximum_accepted_step=std::max(result.maximum_accepted_step,dt);
            save(time,state);
        }else{
            ++result.rejected_steps;
            if(dt<=controls.minimum_step*(1.0+16.0*std::numeric_limits<double>::epsilon()) && remaining>controls.minimum_step)
                throw std::runtime_error("adaptive transient LTE tolerance cannot be met at minimum step");
        }
        const double exponent=use_bdf2?1.0/3.0:0.5;
        const double raw_scale=normalized_error>0.0
            ?controls.safety_factor*std::pow(1.0/normalized_error,exponent)
            :controls.maximum_scale;
        const double factor=std::clamp(raw_scale,controls.minimum_scale,controls.maximum_scale);
        dt=std::clamp(dt*factor,controls.minimum_step,controls.maximum_step);
    }
    if(result.accepted_steps==0U)result.minimum_accepted_step=0.0;
    return result;
}

NoiseResult Circuit::output_noise(double frequency,std::string_view output_node,double temperature,const NewtonConfig& cfg) const {
    if (!(frequency>0.0) || !(temperature>0.0)) throw std::invalid_argument("invalid noise analysis controls");
    const auto operating=dc_operating_point_homotopy(cfg);if(!operating.converged)throw std::runtime_error("noise operating point did not converge");
    const std::size_t output=find_node(output_node),node_unknowns=node_count()-1U,dimension=node_unknowns+branch_count();
    const Complex jw{0.0,2.0*std::numbers::pi*frequency};std::vector<Complex> matrix(dimension*dimension),zero_rhs(dimension);
    for(const auto& r:resistors_) stamp_admittance(matrix,dimension,r.p,r.n,Complex{1.0/r.value,0.0});
    for(const auto& c:capacitors_) stamp_admittance(matrix,dimension,c.p,c.n,jw*c.value);
    for(const auto& s:vccs_) stamp_vccs(matrix,dimension,s.p,s.n,s.cp,s.cn,Complex{s.gain,0});
    for(const auto& source:behavioral_current_sources_){
        const auto linear=behavioral_linearize(source.evaluator,operating.node_voltage,0.0);
        stamp_behavioral_current(matrix,static_cast<std::vector<Complex>*>(nullptr),dimension,source.p,source.n,linear,operating.node_voltage);
    }
    std::size_t branch=0U;
    for(const auto&v:voltage_sources_){
        stamp_voltage_branch(matrix,zero_rhs,dimension,node_count(),branch,v.p,v.n,Complex{});++branch;
    }
    for(const auto& source:behavioral_voltage_sources_){
        const auto linear=behavioral_linearize(source.evaluator,operating.node_voltage,0.0);
        stamp_behavioral_voltage(matrix,static_cast<std::vector<Complex>*>(nullptr),dimension,node_count(),branch,source.p,source.n,linear,operating.node_voltage);
        ++branch;
    }
    for(const auto&e:vcvs_){
        stamp_voltage_branch(matrix,zero_rhs,dimension,node_count(),branch,e.p,e.n,Complex{});
        stamp_vcvs_control(matrix,dimension,node_count(),branch,e.cp,e.cn,Complex{e.gain,0});++branch;
    }
    for(const auto&h:ccvs_){
        stamp_voltage_branch(matrix,zero_rhs,dimension,node_count(),branch,h.p,h.n,Complex{});
        stamp_ccvs_control(matrix,dimension,node_count(),branch,branch_index(h.control),Complex{h.gain,0});++branch;
    }
    for(const auto&l:inductors_){
        stamp_voltage_branch(matrix,zero_rhs,dimension,node_count(),branch,l.p,l.n,Complex{},-jw*l.value);++branch;
    }
    for(const auto&f:cccs_)stamp_cccs(matrix,dimension,node_count(),f.p,f.n,branch_index(f.control),Complex{f.gain,0});
    for(const auto&d:diodes_){const auto&m=diode_models_.at(d.model);const double vt=m.emission_coefficient*boltzmann*m.temperature_k/electron_charge,vd=operating.node_voltage[d.a]-operating.node_voltage[d.k],g=m.saturation_current*std::exp(std::clamp(vd/vt,-40.0,40.0))/vt;stamp_admittance(matrix,dimension,d.a,d.k,Complex{g,0});}
    for(const auto&m:mosfets_){const auto lin=mos_linearize(mos_models_.at(m.model),operating.node_voltage[m.d],operating.node_voltage[m.g],operating.node_voltage[m.s]);stamp_transistor_linearization(matrix,dimension,lin,m.d,m.g,m.s);}for(const auto&j:jfets_){const auto lin=jfet_linearize(jfet_models_.at(j.model),operating.node_voltage[j.d],operating.node_voltage[j.g],operating.node_voltage[j.s]);stamp_transistor_linearization(matrix,dimension,lin,j.d,j.g,j.s);}for(const auto&q:bjts_){const auto lin=bjt_linearize(bjt_models_.at(q.model),operating.node_voltage[q.c],operating.node_voltage[q.b],operating.node_voltage[q.e]);stamp_bjt_jacobian(matrix,dimension,lin,q.c,q.b,q.e);}for(const auto&s:switches_){const double control=operating.node_voltage[s.cp]-operating.node_voltage[s.cn];stamp_admittance(matrix,dimension,s.p,s.n,Complex{switch_conductance(switch_models_.at(s.model),control),0});}
    for(const auto& device:dynamic_devices_)stamp_dynamic_ac(matrix,dimension,device.terminals,device.evaluator,operating.node_voltage,2.0*std::numbers::pi*frequency);
    for(const auto& device:nports_) stamp_nport(matrix,dimension,device,frequency);
    auto transfer=[&](std::size_t p,std::size_t n){std::vector<Complex> rhs(dimension);stamp_current(rhs,p,n,Complex{1.0,0.0});const auto solution=solve_dense(matrix,rhs,dimension);return output?solution[output-1U]:Complex{};};
    NoiseResult result;result.frequency_hz=frequency;
    const auto add_noise=[&](std::string name,std::size_t p,std::size_t n,double source_density){
        if(!(source_density>=0.0)||!std::isfinite(source_density))throw std::runtime_error("non-finite device noise density");
        const Complex h=transfer(p,n);const double contribution=std::norm(h)*source_density;
        result.output_noise_density_v2_per_hz+=contribution;
        result.contributions.push_back({std::move(name),contribution});
    };
    for(const auto&r:resistors_)add_noise(r.name,r.p,r.n,4.0*boltzmann*temperature/r.value);
    for(const auto&d:diodes_){
        const auto&m=diode_models_.at(d.model);const double vt=m.emission_coefficient*boltzmann*m.temperature_k/electron_charge;
        const double vd=operating.node_voltage[d.a]-operating.node_voltage[d.k];
        const double current=m.saturation_current*(std::exp(std::clamp(vd/vt,-40.0,40.0))-1.0);
        add_noise(d.name,d.a,d.k,2.0*electron_charge*std::abs(current));
    }
    for(const auto&m:mosfets_){
        const auto lin=mos_linearize(mos_models_.at(m.model),operating.node_voltage[m.d],operating.node_voltage[m.g],operating.node_voltage[m.s]);
        add_noise(m.name,m.d,m.s,4.0*boltzmann*temperature*(2.0/3.0)*std::max(0.0,lin.gm));
    }
    for(const auto&j:jfets_){
        const auto lin=jfet_linearize(jfet_models_.at(j.model),operating.node_voltage[j.d],operating.node_voltage[j.g],operating.node_voltage[j.s]);
        add_noise(j.name,j.d,j.s,4.0*boltzmann*temperature*(2.0/3.0)*std::max(0.0,lin.gm));
    }
    for(const auto&q:bjts_){
        const auto currents=bjt_currents(bjt_models_.at(q.model),operating.node_voltage[q.c],operating.node_voltage[q.b],operating.node_voltage[q.e]);
        add_noise(q.name+":collector",q.c,q.e,2.0*electron_charge*std::abs(currents[0]));
        add_noise(q.name+":base",q.b,q.e,2.0*electron_charge*std::abs(currents[1]));
    }
    result.output_noise_v_per_sqrt_hz=std::sqrt(std::max(0.0,result.output_noise_density_v2_per_hz));
    return result;
}

double Circuit::voltage(const OperatingPoint& point,std::string_view name) const {const auto id=find_node(name);if(point.node_voltage.size()!=node_count())throw std::invalid_argument("operating point node size mismatch");return point.node_voltage[id];}
Complex Circuit::voltage(const AcSolution& point,std::string_view name) const {const auto id=find_node(name);if(point.node_voltage.size()!=node_count())throw std::invalid_argument("AC point node size mismatch");return point.node_voltage[id];}

Complex Circuit::branch_current(const AcSolution& point,std::string_view name) const {
    if(point.branch_current.size()!=branch_count())throw std::invalid_argument("AC point branch-current size mismatch");
    return point.branch_current[branch_index(name)];
}

double parse_spice_number(std::string_view text) {
    std::string value=upper(trim(std::string(text)));if(value.empty())throw std::invalid_argument("empty SPICE number");std::size_t position=0;double number=std::stod(value,&position);const std::string suffix=value.substr(position);double scale=1.0;
    if(suffix.rfind("MEG",0)==0)scale=1e6;else if(suffix.rfind("T",0)==0)scale=1e12;else if(suffix.rfind("G",0)==0)scale=1e9;else if(suffix.rfind("K",0)==0)scale=1e3;else if(suffix.rfind("M",0)==0)scale=1e-3;else if(suffix.rfind("U",0)==0)scale=1e-6;else if(suffix.rfind("N",0)==0)scale=1e-9;else if(suffix.rfind("P",0)==0)scale=1e-12;else if(suffix.rfind("F",0)==0)scale=1e-15;return number*scale;
}

namespace {
class SpiceExpressionParser {
public:
    SpiceExpressionParser(std::string text,const std::unordered_map<std::string,double>& parameters)
        :text_(std::move(text)),parameters_(parameters) {
        text_=trim(text_);
        if(text_.size()>=2U&&text_.front()=='{'&&text_.back()=='}')text_=text_.substr(1,text_.size()-2U);
    }
    double parse(){position_=0U;const double value=expression();skip();if(position_!=text_.size())throw std::invalid_argument("trailing SPICE expression text");if(!std::isfinite(value))throw std::overflow_error("non-finite SPICE expression");return value;}
private:
    std::string text_;const std::unordered_map<std::string,double>& parameters_;std::size_t position_{};
    void skip(){while(position_<text_.size()&&std::isspace(static_cast<unsigned char>(text_[position_])))++position_;}
    bool take(char c){skip();if(position_<text_.size()&&text_[position_]==c){++position_;return true;}return false;}
    double expression(){double value=term();for(;;){if(take('+'))value+=term();else if(take('-'))value-=term();else return value;}}
    double term(){double value=power();for(;;){if(take('*'))value*=power();else if(take('/')){const double d=power();if(d==0.0)throw std::domain_error("division by zero in SPICE expression");value/=d;}else return value;}}
    double power(){double value=unary();if(take('^'))value=std::pow(value,power());return value;}
    double unary(){if(take('+'))return unary();if(take('-'))return -unary();return primary();}
    double primary(){skip();if(take('(')){const double value=expression();if(!take(')'))throw std::invalid_argument("missing ) in SPICE expression");return value;}if(position_>=text_.size())throw std::invalid_argument("unexpected end of SPICE expression");
        if(std::isalpha(static_cast<unsigned char>(text_[position_]))||text_[position_]=='_'){const std::size_t begin=position_++;while(position_<text_.size()&&(std::isalnum(static_cast<unsigned char>(text_[position_]))||text_[position_]=='_'))++position_;const std::string key=upper(text_.substr(begin,position_-begin));const auto found=parameters_.find(key);if(found==parameters_.end())throw std::out_of_range("unknown SPICE parameter: "+key);return found->second;}
        const std::size_t begin=position_;bool exponent=false;
        while(position_<text_.size()){
            const char c=text_[position_];
            if(std::isdigit(static_cast<unsigned char>(c))||c=='.'){++position_;continue;}
            if((c=='e'||c=='E')&&!exponent){exponent=true;++position_;if(position_<text_.size()&&(text_[position_]=='+'||text_[position_]=='-'))++position_;continue;}
            break;
        }
        if(begin==position_)throw std::invalid_argument("expected SPICE expression value");
        while(position_<text_.size()&&std::isalpha(static_cast<unsigned char>(text_[position_])))++position_;
        return parse_spice_number(text_.substr(begin,position_-begin));
    }
};
}

double evaluate_spice_expression(std::string_view expression,const std::unordered_map<std::string,double>& parameters) {
    std::unordered_map<std::string,double> normalized;normalized.reserve(parameters.size());for(const auto& [name,value]:parameters)normalized[upper(name)]=value;
    return SpiceExpressionParser(std::string(expression),normalized).parse();
}

namespace {
struct SpiceUserFunction {
    std::string name;
    std::vector<std::string> arguments;
    std::string expression;
};

std::string replace_identifier_ci(std::string text,std::string_view identifier,std::string_view replacement) {
    const std::string key=upper(std::string(identifier));
    std::string out;out.reserve(text.size()+replacement.size());
    for(std::size_t i=0;i<text.size();){
        if(std::isalpha(static_cast<unsigned char>(text[i]))||text[i]=='_'){
            const std::size_t begin=i++;
            while(i<text.size()&&(std::isalnum(static_cast<unsigned char>(text[i]))||text[i]=='_'))++i;
            const std::string token=text.substr(begin,i-begin);
            if(upper(token)==key) out.append(replacement); else out.append(token);
        }else out.push_back(text[i++]);
    }
    return out;
}

std::vector<std::string> split_function_arguments(std::string_view text) {
    std::vector<std::string> result;std::size_t begin=0U;int depth=0;
    for(std::size_t i=0;i<=text.size();++i){
        const char c=i<text.size()?text[i]:',';
        if(c=='('||c=='{')++depth; else if(c==')'||c=='}')--depth;
        if(c==','&&depth==0){result.push_back(trim(std::string(text.substr(begin,i-begin))));begin=i+1U;}
    }
    if(result.size()==1U&&result.front().empty())result.clear();
    return result;
}

std::string expand_user_functions(std::string text,
                                  const std::unordered_map<std::string,SpiceUserFunction>& functions,
                                  std::size_t depth=0U) {
    if(depth>32U)throw std::invalid_argument("SPICE .FUNC expansion depth exceeded");
    for(std::size_t i=0;i<text.size();){
        if(!(std::isalpha(static_cast<unsigned char>(text[i]))||text[i]=='_')){++i;continue;}
        const std::size_t begin=i++;
        while(i<text.size()&&(std::isalnum(static_cast<unsigned char>(text[i]))||text[i]=='_'))++i;
        const std::string name=upper(text.substr(begin,i-begin));
        auto found=functions.find(name);
        std::size_t open=i;while(open<text.size()&&std::isspace(static_cast<unsigned char>(text[open])))++open;
        if(found==functions.end()||open>=text.size()||text[open]!='(')continue;
        int nesting=1;std::size_t close=open+1U;
        for(;close<text.size()&&nesting>0;++close){if(text[close]=='(')++nesting;else if(text[close]==')')--nesting;}
        if(nesting!=0)throw std::invalid_argument("unterminated SPICE .FUNC call");
        const std::size_t close_index=close-1U;
        auto actual=split_function_arguments(std::string_view(text).substr(open+1U,close_index-open-1U));
        const auto& function=found->second;
        if(actual.size()!=function.arguments.size())throw std::invalid_argument("SPICE .FUNC argument count mismatch: "+function.name);
        std::string expression=function.expression;
        for(std::size_t a=0;a<actual.size();++a){
            const std::string expanded_actual=expand_user_functions(actual[a],functions,depth+1U);
            expression=replace_identifier_ci(std::move(expression),function.arguments[a],"("+expanded_actual+")");
        }
        expression=expand_user_functions(std::move(expression),functions,depth+1U);
        text.replace(begin,close_index-begin+1U,"("+expression+")");
        i=begin;
    }
    return text;
}

std::vector<std::string> preprocess_user_functions(const std::vector<std::string>& lines) {
    std::unordered_map<std::string,SpiceUserFunction> functions;
    std::vector<std::string> body;body.reserve(lines.size());
    for(const auto& raw:lines){
        const std::string line=trim(raw);auto t=tokens(line);
        if(t.empty()||upper(t[0])!=".FUNC"){body.push_back(raw);continue;}
        const auto name_start=line.find_first_not_of(" \t",5U);
        if(name_start==std::string::npos)throw std::invalid_argument("SPICE .FUNC missing definition");
        const auto open=line.find('(',name_start),close=open==std::string::npos?std::string::npos:line.find(')',open+1U);
        if(open==std::string::npos||close==std::string::npos)throw std::invalid_argument("malformed SPICE .FUNC signature");
        SpiceUserFunction f;f.name=upper(trim(line.substr(name_start,open-name_start)));
        f.arguments=split_function_arguments(std::string_view(line).substr(open+1U,close-open-1U));
        for(auto& arg:f.arguments)arg=upper(arg);
        f.expression=trim(line.substr(close+1U));
        if(f.expression.size()>=2U&&f.expression.front()=='{'&&f.expression.back()=='}')f.expression=f.expression.substr(1U,f.expression.size()-2U);
        if(f.name.empty()||f.expression.empty())throw std::invalid_argument("malformed SPICE .FUNC definition");
        functions[f.name]=std::move(f);
    }
    if(functions.empty())return body;
    for(auto& line:body)line=expand_user_functions(std::move(line),functions);
    return body;
}

std::string unquote_path(std::string token) {
    token=trim(std::move(token));
    if(token.size()>=2U&&((token.front()=='\"'&&token.back()=='\"')||(token.front()=='\''&&token.back()=='\'')))
        return token.substr(1U,token.size()-2U);
    return token;
}

std::vector<std::string> load_spice_file_recursive(const std::filesystem::path& input_path,
                                                   std::string_view requested_section,
                                                   std::unordered_set<std::string>& active) {
    std::error_code ec;const auto canonical=std::filesystem::weakly_canonical(input_path,ec);
    const auto path=ec?input_path:canonical;const std::string key=path.string();
    if(!active.insert(key).second)throw std::invalid_argument("recursive SPICE include: "+key);
    std::ifstream stream(path);if(!stream)throw std::runtime_error("cannot open SPICE file: "+key);
    std::vector<std::string> raw;std::string line;while(std::getline(stream,line))raw.push_back(line);
    const auto base=path.parent_path();std::vector<std::string> selected;
    bool section_active=requested_section.empty(),found_section=requested_section.empty();int library_depth=0;
    for(const auto& raw_line:raw){
        const std::string stripped=trim(raw_line);auto t=tokens(stripped);
        if(t.empty()){if(section_active)selected.push_back(raw_line);continue;}
        const std::string directive=upper(t[0]);
        if(!requested_section.empty()){
            if(directive==".LIB"&&t.size()==2U&&upper(t[1])==upper(std::string(requested_section))){section_active=true;found_section=true;continue;}
            if(section_active&&directive==".ENDL"){section_active=false;continue;}
            if(!section_active)continue;
        }else{
            if(directive==".LIB"&&t.size()==2U){++library_depth;continue;}
            if(directive==".ENDL"&&library_depth>0){--library_depth;continue;}
            if(library_depth>0)continue;
        }
        if(directive==".INCLUDE"){
            if(t.size()<2U)throw std::invalid_argument("SPICE .INCLUDE requires a path");
            auto nested=load_spice_file_recursive(base/unquote_path(t[1]),{},active);selected.insert(selected.end(),nested.begin(),nested.end());continue;
        }
        if(directive==".LIB"&&t.size()>=3U){
            auto nested=load_spice_file_recursive(base/unquote_path(t[1]),t[2],active);selected.insert(selected.end(),nested.begin(),nested.end());continue;
        }
        selected.push_back(raw_line);
    }
    active.erase(key);
    if(!found_section)throw std::invalid_argument("SPICE .LIB section not found: "+std::string(requested_section));
    return selected;
}

struct SubcircuitDefinition {
    std::string name;
    std::vector<std::string> pins;
    std::unordered_map<std::string,std::string> default_parameters;
    std::vector<std::string> body;
};

std::string numeric_text(double value) {
    std::ostringstream out;
    out<<std::setprecision(17)<<value;
    return out.str();
}

std::unordered_map<std::string,double> parse_parameter_assignments(
    std::span<const std::string> assignments,
    const std::unordered_map<std::string,double>& inherited) {
    auto parameters=inherited;
    for (const auto& token:assignments) {
        if (upper(token)=="PARAMS:") continue;
        const auto equal=token.find('=');
        if (equal==std::string::npos) continue;
        const std::string key=upper(token.substr(0,equal));
        parameters[key]=evaluate_spice_expression(token.substr(equal+1U),parameters);
    }
    return parameters;
}

std::string map_subcircuit_node(std::string token,
                                const std::unordered_map<std::string,std::string>& node_map,
                                std::string_view prefix) {
    if (token=="0") return token;
    const auto found=node_map.find(upper(token));
    if (found!=node_map.end()) return found->second;
    return std::string(prefix)+token;
}

std::string substitute_parameter_token(std::string token,
                                       const std::unordered_map<std::string,double>& parameters) {
    if (token.empty()) return token;
    const auto exact=parameters.find(upper(token));
    if (exact!=parameters.end()) return numeric_text(exact->second);
    std::size_t begin=0U;
    while ((begin=token.find('{',begin))!=std::string::npos) {
        const auto end=token.find('}',begin+1U);
        if (end==std::string::npos) throw std::invalid_argument("unclosed parameter expression in subcircuit");
        const double value=evaluate_spice_expression(token.substr(begin,end-begin+1U),parameters);
        token.replace(begin,end-begin+1U,numeric_text(value));
        begin+=1U;
    }
    return token;
}

std::string map_behavior_voltage_references(std::string text,
                                            const std::unordered_map<std::string,std::string>& node_map,
                                            std::string_view prefix) {
    std::string out;out.reserve(text.size()+16U);
    for(std::size_t i=0;i<text.size();){
        if((text[i]=='V'||text[i]=='v')&&i+1U<text.size()&&text[i+1U]=='('){
            const auto close=text.find(')',i+2U);
            if(close==std::string::npos)throw std::invalid_argument("unclosed V() reference in behavioral source");
            const auto args=split_function_arguments(std::string_view(text).substr(i+2U,close-i-2U));
            if(args.empty()||args.size()>2U)throw std::invalid_argument("behavioral V() requires one or two node names");
            out+="V(";
            out+=map_subcircuit_node(args[0],node_map,prefix);
            if(args.size()==2U){out+=',';out+=map_subcircuit_node(args[1],node_map,prefix);}
            out+=')';i=close+1U;continue;
        }
        out.push_back(text[i++]);
    }
    return out;
}

std::string substitute_behavior_parameters(std::string text,
                                           const std::unordered_map<std::string,double>& parameters) {
    std::string out;out.reserve(text.size()+16U);
    for(std::size_t i=0;i<text.size();){
        if((text[i]=='V'||text[i]=='v')&&i+1U<text.size()&&text[i+1U]=='('){
            const auto close=text.find(')',i+2U);
            if(close==std::string::npos)throw std::invalid_argument("unclosed V() reference in behavioral source");
            out.append(text,i,close-i+1U);i=close+1U;continue;
        }
        if(std::isalpha(static_cast<unsigned char>(text[i]))||text[i]=='_'){
            const std::size_t begin=i++;
            while(i<text.size()&&(std::isalnum(static_cast<unsigned char>(text[i]))||text[i]=='_'))++i;
            const std::string token=text.substr(begin,i-begin);
            const auto found=parameters.find(upper(token));
            if(found!=parameters.end())out+=numeric_text(found->second);else out+=token;
            continue;
        }
        out.push_back(text[i++]);
    }
    return out;
}

std::string flatten_device_line(const std::string& raw,
                                const std::unordered_map<std::string,std::string>& node_map,
                                const std::unordered_map<std::string,double>& parameters,
                                const std::unordered_map<std::string,std::string>& model_map,
                                std::string_view prefix) {
    auto t=tokens(raw);
    if (t.empty()) return {};
    const char kind=static_cast<char>(std::toupper(static_cast<unsigned char>(t[0][0])));
    const auto prefixed_element=[&](const std::string& name){
        if(name.empty()||prefix.empty()) return name;
        return name.substr(0,1)+std::string(prefix)+name.substr(1);
    };
    t[0]=prefixed_element(t[0]);
    auto map_nodes=[&](std::size_t count){
        for(std::size_t i=1U;i<=count&&i<t.size();++i)t[i]=map_subcircuit_node(t[i],node_map,prefix);
    };
    const std::string dependent_mode=(kind=='G'||kind=='E')&&t.size()>3U?upper(t[3]):std::string{};
    const bool dependent_poly=dependent_mode.rfind("POLY(",0U)==0U;
    const bool dependent_table=dependent_mode=="TABLE";
    if(kind=='R'||kind=='C'||kind=='L'||kind=='D'||kind=='V'||kind=='I'||kind=='B') map_nodes(2U);
    else if((kind=='G'||kind=='E')&&(dependent_poly||dependent_table)){
        map_nodes(2U);
        if(dependent_poly&&t.size()>=6U){
            t[4]=map_subcircuit_node(t[4],node_map,prefix);
            t[5]=map_subcircuit_node(t[5],node_map,prefix);
        }
    }
    else if(kind=='G'||kind=='E'||kind=='S') map_nodes(4U);
    else if(kind=='F'||kind=='H') map_nodes(2U);
    else if(kind=='M') map_nodes(4U);
    else if(kind=='Q'||kind=='J') map_nodes(3U);
    if(kind=='K'&&t.size()>=4U){t[1]=prefixed_element(t[1]);t[2]=prefixed_element(t[2]);}
    if((kind=='F'||kind=='H')&&t.size()>=4U)t[3]=prefixed_element(t[3]);
    std::size_t model_index=static_cast<std::size_t>(-1);
    if(kind=='D')model_index=3U;else if(kind=='M')model_index=5U;else if(kind=='Q'||kind=='J')model_index=4U;else if(kind=='S')model_index=5U;
    if(model_index<t.size()){const auto model=model_map.find(upper(t[model_index]));if(model!=model_map.end())t[model_index]=model->second;}

    // Only substitute parameter-bearing value fields; model and node names stay literal.
    if((kind=='R'||kind=='C'||kind=='L')&&t.size()>=4U)t[3]=substitute_parameter_token(t[3],parameters);
    else if((kind=='G'||kind=='E')&&dependent_poly){
        for(std::size_t i=6U;i<t.size();++i)t[i]=substitute_parameter_token(t[i],parameters);
    }
    else if((kind=='G'||kind=='E')&&dependent_table){
        std::string specification;for(std::size_t i=4U;i<t.size();++i){if(i>4U)specification+=' ';specification+=t[i];}
        specification=map_behavior_voltage_references(std::move(specification),node_map,prefix);
        specification=substitute_behavior_parameters(std::move(specification),parameters);
        t.resize(5U);t[4]=std::move(specification);
    }
    else if((kind=='G'||kind=='E')&&t.size()>=6U)t[5]=substitute_parameter_token(t[5],parameters);
    else if((kind=='F'||kind=='H')&&t.size()>=5U)t[4]=substitute_parameter_token(t[4],parameters);
    else if(kind=='K'&&t.size()>=4U)t[3]=substitute_parameter_token(t[3],parameters);
    else if(kind=='V'||kind=='I')for(std::size_t i=3U;i<t.size();++i)t[i]=substitute_parameter_token(t[i],parameters);
    else if(kind=='B'&&t.size()>=4U){
        std::string specification;
        for(std::size_t i=3U;i<t.size();++i){if(i>3U)specification+=' ';specification+=t[i];}
        specification=map_behavior_voltage_references(std::move(specification),node_map,prefix);
        specification=substitute_behavior_parameters(std::move(specification),parameters);
        t.resize(4U);t[3]=std::move(specification);
    }

    std::ostringstream out;
    for(std::size_t i=0;i<t.size();++i){if(i)out<<' ';out<<t[i];}
    return out.str();
}

void expand_subcircuit_instance(const std::vector<std::string>& instance,
                                const std::unordered_map<std::string,SubcircuitDefinition>& definitions,
                                const std::unordered_map<std::string,double>& inherited_parameters,
                                const std::unordered_map<std::string,std::string>& outer_nodes,
                                const std::unordered_map<std::string,std::string>& inherited_models,
                                std::string_view outer_prefix,
                                std::vector<std::string>& output,
                                std::vector<std::string>& stack) {
    const SubcircuitDefinition* definition=nullptr;
    std::size_t name_index=0U;
    for(const auto& [name,candidate]:definitions){
        const std::size_t index=1U+candidate.pins.size();
        if(index<instance.size()&&upper(instance[index])==name){definition=&candidate;name_index=index;break;}
    }
    if(!definition)throw std::invalid_argument("unknown or malformed SPICE subcircuit instance: "+instance.front());
    if(std::find(stack.begin(),stack.end(),definition->name)!=stack.end())throw std::invalid_argument("recursive SPICE subcircuit hierarchy");

    std::unordered_map<std::string,double> parameters=inherited_parameters;
    for(const auto& [key,expression]:definition->default_parameters)
        parameters[key]=evaluate_spice_expression(expression,parameters);
    parameters=parse_parameter_assignments(std::span<const std::string>(instance).subspan(name_index+1U),parameters);

    std::unordered_map<std::string,std::string> node_map;
    for(std::size_t i=0;i<definition->pins.size();++i)
        node_map[definition->pins[i]]=map_subcircuit_node(instance[i+1U],outer_nodes,outer_prefix);
    const std::string prefix=std::string(outer_prefix)+instance.front()+":";
    auto models=inherited_models;
    for(const auto& body_line:definition->body){auto mt=tokens(body_line);if(mt.size()>=3U&&upper(mt[0])==".MODEL")models[upper(mt[1])]=prefix+mt[1];}
    stack.push_back(definition->name);
    for(const auto& body_line:definition->body){
        auto t=tokens(body_line);if(t.empty())continue;
        if(t[0][0]=='.'&&upper(t[0])==".PARAM"){
            parameters=parse_parameter_assignments(std::span<const std::string>(t).subspan(1U),parameters);
            continue;
        }
        if(t[0][0]=='.'&&upper(t[0])==".MODEL"){
            if(t.size()<3U)throw std::invalid_argument("malformed local SPICE .MODEL");
            t[1]=models.at(upper(t[1]));
            for(std::size_t k=2U;k<t.size();++k)t[k]=substitute_parameter_token(t[k],parameters);
            std::ostringstream model_line;for(std::size_t k=0;k<t.size();++k){if(k)model_line<<' ';model_line<<t[k];}output.push_back(model_line.str());
            continue;
        }
        if(std::toupper(static_cast<unsigned char>(t[0][0]))=='X'){
            expand_subcircuit_instance(t,definitions,parameters,node_map,models,prefix,output,stack);
        }else output.push_back(flatten_device_line(body_line,node_map,parameters,models,prefix));
    }
    stack.pop_back();
}

std::vector<std::string> flatten_subcircuits(const std::vector<std::string>& lines) {
    std::unordered_map<std::string,SubcircuitDefinition> definitions;
    std::vector<std::string> top_level;
    for(std::size_t i=0;i<lines.size();++i){
        std::string line=trim(lines[i]);
        if(line.empty()||line[0]=='*')continue;
        auto t=tokens(line);
        if(!t.empty()&&upper(t[0])==".SUBCKT"){
            if(t.size()<3U)throw std::invalid_argument("SPICE .SUBCKT requires name and pins");
            SubcircuitDefinition definition;definition.name=upper(t[1]);
            std::size_t parameter_start=t.size();
            for(std::size_t k=2U;k<t.size();++k)if(upper(t[k])=="PARAMS:"||t[k].find('=')!=std::string::npos){parameter_start=k;break;}
            for(std::size_t k=2U;k<parameter_start;++k)definition.pins.push_back(upper(t[k]));
            for(std::size_t k=parameter_start;k<t.size();++k){
                if(upper(t[k])=="PARAMS:") continue;
                const auto equal=t[k].find('=');
                if(equal!=std::string::npos) definition.default_parameters[upper(t[k].substr(0,equal))]=t[k].substr(equal+1U);
            }
            bool ended=false;
            for(++i;i<lines.size();++i){
                std::string body=trim(lines[i]);auto bt=tokens(body);
                if(!bt.empty()&&upper(bt[0])==".ENDS"){ended=true;break;}
                definition.body.push_back(body);
            }
            if(!ended)throw std::invalid_argument("unterminated SPICE .SUBCKT "+definition.name);
            definitions[definition.name]=std::move(definition);
        }else top_level.push_back(line);
    }
    if(definitions.empty())return top_level;

    std::unordered_map<std::string,double> global_parameters;
    for(const auto& line:top_level){auto t=tokens(line);if(!t.empty()&&upper(t[0])==".PARAM")global_parameters=parse_parameter_assignments(std::span<const std::string>(t).subspan(1U),global_parameters);}
    std::vector<std::string> flattened;std::vector<std::string> stack;
    const std::unordered_map<std::string,std::string> no_nodes,no_models;
    for(const auto& line:top_level){
        auto t=tokens(line);if(t.empty())continue;
        if(std::toupper(static_cast<unsigned char>(t[0][0]))=='X')expand_subcircuit_instance(t,definitions,global_parameters,no_nodes,no_models,"",flattened,stack);
        else flattened.push_back(line);
    }
    return flattened;
}
} // namespace

Circuit Circuit::parse_spice(std::istream& input) {
    Circuit circuit;
    std::string input_line;
    std::vector<std::string> raw_lines;
    while(std::getline(input,input_line)) raw_lines.push_back(input_line);
    for(const auto& raw:raw_lines){auto t=tokens(trim(raw));if(!t.empty()&&(upper(t[0])==".INCLUDE"||(upper(t[0])==".LIB"&&t.size()>=3U)))throw std::invalid_argument("file directives require Circuit::parse_spice_file");}
    const auto functional_lines=preprocess_user_functions(raw_lines);
    const auto flattened_lines=flatten_subcircuits(functional_lines);
    std::vector<std::string> devices,mutuals,models,param_lines,ic_lines,nodeset_lines;
    for(std::string line:flattened_lines){
        line=trim(line);
        if(line.empty()||line[0]=='*')continue;
        if(line[0]=='.'){
            const auto t=tokens(line);
            if(!t.empty()&&upper(t[0])==".MODEL")models.push_back(line);
            else if(!t.empty()&&upper(t[0])==".PARAM")param_lines.push_back(line);
            else if(!t.empty()&&upper(t[0])==".IC")ic_lines.push_back(line);
            else if(!t.empty()&&upper(t[0])==".NODESET")nodeset_lines.push_back(line);
            continue;
        }
        if(std::toupper(static_cast<unsigned char>(line[0]))=='K')mutuals.push_back(line);
        else devices.push_back(line);
    }

    std::unordered_map<std::string,double> parameters;
    for(auto raw:param_lines){
        for(char& c:raw)if(c==',')c=' ';
        const auto t=tokens(raw);
        for(std::size_t i=1;i<t.size();++i){
            const auto equal=t[i].find('=');
            if(equal==std::string::npos)continue;
            const std::string name=upper(t[i].substr(0,equal));
            parameters[name]=evaluate_spice_expression(t[i].substr(equal+1U),parameters);
        }
    }
    const auto value=[&](std::string_view token){return evaluate_spice_expression(token,parameters);};
    const auto behavioral_evaluator=[&](std::string expression)->BehavioralEvaluator{
        expression=trim(std::move(expression));
        if(expression.empty())throw std::invalid_argument("empty behavioral-source expression");
        struct VoltageReference { std::size_t positive{},negative{}; std::string parameter; };
        std::vector<VoltageReference> references;
        std::string transformed;transformed.reserve(expression.size()+16U);
        for(std::size_t i=0;i<expression.size();){
            if((expression[i]=='V'||expression[i]=='v')&&i+1U<expression.size()&&expression[i+1U]=='('){
                const auto close=expression.find(')',i+2U);
                if(close==std::string::npos)throw std::invalid_argument("unclosed V() reference in behavioral source");
                const auto args=split_function_arguments(std::string_view(expression).substr(i+2U,close-i-2U));
                if(args.empty()||args.size()>2U)throw std::invalid_argument("behavioral V() requires one or two node names");
                VoltageReference reference;
                reference.positive=circuit.node(args[0]);
                reference.negative=args.size()==2U?circuit.node(args[1]):0U;
                reference.parameter="__BV"+std::to_string(references.size());
                transformed+=reference.parameter;references.push_back(std::move(reference));i=close+1U;continue;
            }
            transformed.push_back(expression[i++]);
        }
        const auto captured_parameters=parameters;
        return [transformed=std::move(transformed),references=std::move(references),captured_parameters]
            (std::span<const double> node_voltage,double time_s){
                std::unordered_map<std::string,double> local=captured_parameters;
                local["TIME"]=time_s;
                for(const auto& reference:references){
                    if(reference.positive>=node_voltage.size()||reference.negative>=node_voltage.size())
                        throw std::out_of_range("behavioral-source node reference out of range");
                    local[reference.parameter]=node_voltage[reference.positive]-node_voltage[reference.negative];
                }
                return evaluate_spice_expression(transformed,local);
            };
    };

    const auto source_waveform=[&](std::string specification)->std::pair<double,SourceWaveform>{
        for(char& c:specification) if(c=='('||c==')'||c==',') c=' ';
        const auto args=tokens(specification);
        if(args.empty()) throw std::invalid_argument("empty transient source specification");
        const std::string kind=upper(args.front());
        if(kind=="PULSE") {
            if(args.size()<8U) throw std::invalid_argument("PULSE requires V1 V2 TD TR TF PW PER");
            PulseWaveform pulse{value(args[1]),value(args[2]),value(args[3]),value(args[4]),
                                value(args[5]),value(args[6]),value(args[7])};
            return {pulse.initial,[pulse](double time){return pulse(time);}};
        }
        if(kind=="SIN") {
            if(args.size()<4U) throw std::invalid_argument("SIN requires VO VA FREQ");
            SineWaveform sine;
            sine.offset=value(args[1]);sine.amplitude=value(args[2]);sine.frequency_hz=value(args[3]);
            if(args.size()>4U && upper(args[4])!="AC") sine.delay=value(args[4]);
            if(args.size()>5U && upper(args[5])!="AC") sine.damping_per_s=value(args[5]);
            if(args.size()>6U && upper(args[6])!="AC") sine.phase_deg=value(args[6]);
            return {sine.offset,[sine](double time){return sine(time);}};
        }
        if(kind=="PWL") {
            PwlWaveform pwl;
            for(std::size_t i=1U;i+1U<args.size() && upper(args[i])!="AC";i+=2U)
                pwl.points.emplace_back(value(args[i]),value(args[i+1U]));
            if(pwl.points.empty()) throw std::invalid_argument("PWL requires time/value pairs");
            const double dc=pwl.points.front().second;
            return {dc,[pwl](double time){return pwl(time);}};
        }
        throw std::invalid_argument("unsupported transient source waveform");
    };

    for(const auto& raw:models){
        std::string normalized_model=raw;for(char& c:normalized_model)if(c=='('||c==')'||c==',')c=' ';
        const auto t=tokens(normalized_model);if(t.size()<3U)continue;
        const std::string name=t[1],type=upper(t[2]);const auto params=model_params(raw);
        const auto model_value=[&](std::string_view key,double fallback){const auto found=params.find(std::string(key));return found==params.end()?fallback:value(found->second);};
        if(type=="D"){
            DiodeModel model;model.saturation_current=model_value("IS",model.saturation_current);model.emission_coefficient=model_value("N",model.emission_coefficient);circuit.add_diode_model(name,model);
        }else if(type=="NMOS"||type=="PMOS"){
            MosLevel1Model model;model.pmos=type=="PMOS";model.threshold_voltage=model_value("VTO",model.threshold_voltage);model.transconductance_parameter=model_value("KP",model.transconductance_parameter);model.channel_length_modulation=model_value("LAMBDA",model.channel_length_modulation);circuit.add_mos_model(name,model);
        }else if(type=="NPN"||type=="PNP"){
            BjtModel model;model.pnp=type=="PNP";model.saturation_current=model_value("IS",model.saturation_current);model.forward_beta=model_value("BF",model.forward_beta);model.reverse_beta=model_value("BR",model.reverse_beta);circuit.add_bjt_model(name,model);
        }else if(type=="NJF"||type=="PJF"){
            JfetModel model;model.pchannel=type=="PJF";model.pinch_off_voltage=model_value("VTO",model.pinch_off_voltage);model.beta=model_value("BETA",model.beta);model.channel_length_modulation=model_value("LAMBDA",model.channel_length_modulation);circuit.add_jfet_model(name,model);
        }else if(type=="SW"){
            SwitchModel model;model.on_resistance=model_value("RON",model.on_resistance);model.off_resistance=model_value("ROFF",model.off_resistance);model.threshold_voltage=model_value("VT",model.threshold_voltage);model.hysteresis_voltage=model_value("VH",model.hysteresis_voltage);circuit.add_switch_model(name,model);
        }
    }

    for(const auto& raw:devices){
        const auto t=tokens(raw);if(t.empty())continue;
        const char kind=static_cast<char>(std::toupper(static_cast<unsigned char>(t[0][0])));
        if((kind=='R'||kind=='C'||kind=='L')&&t.size()>=4U){
            const auto p=circuit.node(t[1]),n=circuit.node(t[2]);const double numeric=value(t[3]);
            if(kind=='R')circuit.add_resistor(t[0],p,n,numeric);else if(kind=='C')circuit.add_capacitor(t[0],p,n,numeric);else circuit.add_inductor(t[0],p,n,numeric);
        }else if((kind=='V'||kind=='I')&&t.size()>=4U){
            const auto p=circuit.node(t[1]),n=circuit.node(t[2]);
            double dc=0.0;Complex ac_value{};SourceWaveform transient_waveform;
            std::string specification;
            for(std::size_t k=3U;k<t.size();++k){if(k>3U)specification+=' ';specification+=t[k];}
            const std::string first=upper(t[3]);
            if(first.rfind("PULSE",0)==0U || first.rfind("SIN",0)==0U || first.rfind("PWL",0)==0U){
                const auto parsed=source_waveform(specification);dc=parsed.first;transient_waveform=parsed.second;
            } else {
                std::size_t i=3U;
                if(upper(t[i])=="DC"&&i+1U<t.size()){dc=value(t[++i]);++i;}else{dc=value(t[i]);++i;}
            }
            for(std::size_t i=3U;i<t.size();++i)if(upper(t[i])=="AC"&&i+1U<t.size()){
                const double magnitude=value(t[++i]);double phase=0.0;
                if(i+1U<t.size())phase=value(t[++i]);
                ac_value=std::polar(magnitude,phase*std::numbers::pi/180.0);break;
            }
            if(kind=='V')circuit.add_voltage_source(t[0],p,n,dc,ac_value,std::move(transient_waveform));
            else circuit.add_current_source(t[0],p,n,dc,ac_value,std::move(transient_waveform));
        }else if(kind=='B'&&t.size()>=4U){
            const auto p=circuit.node(t[1]),n=circuit.node(t[2]);
            std::string specification;for(std::size_t i=3U;i<t.size();++i){if(i>3U)specification+=' ';specification+=t[i];}
            const auto equal=specification.find('=');
            if(equal==std::string::npos)throw std::invalid_argument("behavioral source requires V={expr} or I={expr}");
            const std::string output=upper(trim(specification.substr(0U,equal)));
            auto evaluator=behavioral_evaluator(specification.substr(equal+1U));
            if(output=="V")circuit.add_behavioral_voltage_source(t[0],p,n,std::move(evaluator));
            else if(output=="I")circuit.add_behavioral_current_source(t[0],p,n,std::move(evaluator));
            else throw std::invalid_argument("behavioral source output must be V or I");
        }else if((kind=='E'||kind=='G')&&t.size()>=4U){
            const auto p=circuit.node(t[1]),n=circuit.node(t[2]);
            const std::string mode=upper(t[3]);
            if(mode.rfind("POLY(",0U)==0U){
                const auto close=mode.find(')');
                if(close==std::string::npos)throw std::invalid_argument("malformed dependent-source POLY() declaration");
                const auto dimensions=static_cast<std::size_t>(std::stoul(mode.substr(5U,close-5U)));
                if(dimensions!=1U)throw std::invalid_argument("dependent-source POLY baseline currently supports POLY(1)");
                if(t.size()<7U)throw std::invalid_argument("POLY(1) source requires control pair and coefficients");
                const auto cp=circuit.node(t[4]),cn=circuit.node(t[5]);
                std::vector<double> coefficients;coefficients.reserve(t.size()-6U);
                for(std::size_t i=6U;i<t.size();++i)coefficients.push_back(value(t[i]));
                BehavioralEvaluator evaluator=[cp,cn,coefficients=std::move(coefficients)](std::span<const double> nodes,double){
                    const double x=nodes[cp]-nodes[cn];double result=0.0;
                    for(auto it=coefficients.rbegin();it!=coefficients.rend();++it)result=result*x+*it;
                    return result;
                };
                if(kind=='E')circuit.add_behavioral_voltage_source(t[0],p,n,std::move(evaluator));
                else circuit.add_behavioral_current_source(t[0],p,n,std::move(evaluator));
            }else if(mode=="TABLE"){
                std::string specification;for(std::size_t i=4U;i<t.size();++i){if(i>4U)specification+=' ';specification+=t[i];}
                const auto equal=specification.find('=');
                if(equal==std::string::npos)throw std::invalid_argument("TABLE source requires input-expression = (x,y) points");
                auto table_input=behavioral_evaluator(trim(specification.substr(0U,equal)));
                std::string point_text=specification.substr(equal+1U);
                for(char& c:point_text)if(c=='('||c==')'||c==',')c=' ';
                const auto point_tokens=tokens(point_text);
                if(point_tokens.size()<4U||point_tokens.size()%2U!=0U)
                    throw std::invalid_argument("TABLE source requires at least two (x,y) points");
                std::vector<std::pair<double,double>> points;points.reserve(point_tokens.size()/2U);
                for(std::size_t i=0;i<point_tokens.size();i+=2U)points.emplace_back(value(point_tokens[i]),value(point_tokens[i+1U]));
                std::sort(points.begin(),points.end(),[](const auto& a,const auto& b){return a.first<b.first;});
                for(std::size_t i=1;i<points.size();++i)if(!(points[i].first>points[i-1U].first))
                    throw std::invalid_argument("TABLE abscissas must be strictly increasing");
                BehavioralEvaluator evaluator=[table_input=std::move(table_input),points=std::move(points)](std::span<const double> nodes,double time){
                    const double x=table_input(nodes,time);
                    if(x<=points.front().first)return points.front().second;
                    if(x>=points.back().first)return points.back().second;
                    const auto upper_point=std::upper_bound(points.begin(),points.end(),x,
                        [](double sample,const auto& point){return sample<point.first;});
                    const auto lower_point=std::prev(upper_point);
                    const double fraction=(x-lower_point->first)/(upper_point->first-lower_point->first);
                    return lower_point->second+fraction*(upper_point->second-lower_point->second);
                };
                if(kind=='E')circuit.add_behavioral_voltage_source(t[0],p,n,std::move(evaluator));
                else circuit.add_behavioral_current_source(t[0],p,n,std::move(evaluator));
            }else{
                if(t.size()<6U)throw std::invalid_argument("dependent source requires control nodes and gain");
                if(kind=='G')circuit.add_vccs(t[0],p,n,circuit.node(t[3]),circuit.node(t[4]),value(t[5]));
                else circuit.add_vcvs(t[0],p,n,circuit.node(t[3]),circuit.node(t[4]),value(t[5]));
            }
        }else if(kind=='F'&&t.size()>=5U)circuit.add_cccs(t[0],circuit.node(t[1]),circuit.node(t[2]),t[3],value(t[4]));
        else if(kind=='H'&&t.size()>=5U)circuit.add_ccvs(t[0],circuit.node(t[1]),circuit.node(t[2]),t[3],value(t[4]));
        else if(kind=='D'&&t.size()>=4U)circuit.add_diode(t[0],circuit.node(t[1]),circuit.node(t[2]),t[3]);
        else if(kind=='M'&&t.size()>=6U)circuit.add_mosfet(t[0],circuit.node(t[1]),circuit.node(t[2]),circuit.node(t[3]),circuit.node(t[4]),t[5]);
        else if(kind=='Q'&&t.size()>=5U)circuit.add_bjt(t[0],circuit.node(t[1]),circuit.node(t[2]),circuit.node(t[3]),t.back());
        else if(kind=='J'&&t.size()>=5U)circuit.add_jfet(t[0],circuit.node(t[1]),circuit.node(t[2]),circuit.node(t[3]),t[4]);
        else if(kind=='S'&&t.size()>=6U)circuit.add_switch(t[0],circuit.node(t[1]),circuit.node(t[2]),circuit.node(t[3]),circuit.node(t[4]),t[5]);
    }
    for(const auto& raw:mutuals){const auto t=tokens(raw);if(t.size()>=4U)circuit.add_mutual_inductance(t[0],t[1],t[2],value(t[3]));}
    const auto apply_voltage_directive=[&](const std::string& raw,bool initial){
        const auto first_space=raw.find_first_of(" \t");
        if(first_space==std::string::npos)throw std::invalid_argument("empty SPICE voltage directive");
        const std::string body=raw.substr(first_space+1U);
        std::size_t position=0U;bool any=false;
        while(position<body.size()){
            while(position<body.size()&&std::isspace(static_cast<unsigned char>(body[position])))++position;
            if(position>=body.size())break;
            if(std::toupper(static_cast<unsigned char>(body[position]))!='V' || position+1U>=body.size() || body[position+1U]!='(')
                throw std::invalid_argument("SPICE .IC/.NODESET currently expects V(node)=value assignments");
            const auto close=body.find(')',position+2U);
            if(close==std::string::npos)throw std::invalid_argument("malformed SPICE voltage directive node");
            const std::string node_name=trim(body.substr(position+2U,close-position-2U));
            if(node_name.empty()||node_name.find(',')!=std::string::npos)
                throw std::invalid_argument("SPICE .IC/.NODESET differential assignments are not supported yet");
            std::size_t equal=close+1U;while(equal<body.size()&&std::isspace(static_cast<unsigned char>(body[equal])))++equal;
            if(equal>=body.size()||body[equal]!='=')throw std::invalid_argument("SPICE voltage directive requires =");
            std::size_t begin=equal+1U;while(begin<body.size()&&std::isspace(static_cast<unsigned char>(body[begin])))++begin;
            std::size_t end=begin;int depth=0;
            for(;end<body.size();++end){
                const char c=body[end];if(c=='{'||c=='(')++depth;else if(c=='}'||c==')')--depth;
                if(depth==0&&std::isspace(static_cast<unsigned char>(c)))break;
            }
            if(begin==end)throw std::invalid_argument("SPICE voltage directive missing value");
            const double voltage=value(body.substr(begin,end-begin));
            const auto node_id=circuit.node(node_name);
            if(initial)circuit.set_initial_voltage(node_id,voltage);else circuit.set_nodeset_voltage(node_id,voltage);
            any=true;position=end;
        }
        if(!any)throw std::invalid_argument("empty SPICE voltage directive");
    };
    for(const auto& raw:ic_lines)apply_voltage_directive(raw,true);
    for(const auto& raw:nodeset_lines)apply_voltage_directive(raw,false);
    return circuit;
}

Circuit Circuit::parse_spice(std::string_view text){std::istringstream input{std::string(text)};return parse_spice(input);}

Circuit Circuit::parse_spice_file(const std::string& path,std::string_view library_section){
    std::unordered_set<std::string> active;
    const auto lines=load_spice_file_recursive(std::filesystem::path(path),library_section,active);
    std::ostringstream merged;for(const auto& line:lines)merged<<line<<'\n';
    return parse_spice(merged.str());
}

} // namespace cfd::circuit

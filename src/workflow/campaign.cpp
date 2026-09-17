#include "cfd/workflow/campaign.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace cfd::workflow {
namespace {
std::string case_id(std::size_t index){std::ostringstream out;out<<"case"<<std::setw(4)<<std::setfill('0')<<(index+1U);return out.str();}
bool truthy(std::string_view value){return !(value.empty()||value=="false"||value=="False"||value=="FALSE"||value=="0");}
std::string number_to_string(double value){std::ostringstream out;out<<std::setprecision(12)<<value;return out.str();}
void validate_parameter(const CampaignParameter& p){
    if(p.name.empty())throw std::invalid_argument("campaign parameter requires a name");
    if(p.name=="case_id")throw std::invalid_argument("case_id is reserved");
    if(p.kind==CampaignParameterKind::continuous){if(!(p.maximum>=p.minimum)||!std::isfinite(p.minimum)||!std::isfinite(p.maximum))throw std::invalid_argument("invalid continuous campaign parameter");}
    else if(p.levels.empty())throw std::invalid_argument("discrete campaign parameter requires levels");
}
}

double CampaignRandom::uniform_open01(){state=state*2862933555777941757ULL+3037000493ULL;const std::uint64_t mantissa=(state>>11)|1ULL;return std::min(std::nextafter(1.0,0.0),std::max(std::numeric_limits<double>::min(),static_cast<double>(mantissa)*(1.0/9007199254740992.0)));}

std::vector<CampaignCase> generate_factorial_campaign(std::span<const CampaignParameter> parameters){
    if(parameters.empty())throw std::invalid_argument("factorial campaign requires parameters");
    for(const auto& p:parameters){validate_parameter(p);if(p.kind!=CampaignParameterKind::discrete)throw std::invalid_argument("factorial campaign only supports discrete parameters");}
    std::vector<CampaignCase> cases(1U);cases[0].case_id=case_id(0U);
    for(const auto& parameter:parameters){
        std::vector<CampaignCase> next;next.reserve(cases.size()*parameter.levels.size());
        for(const auto& base:cases)for(const auto& level:parameter.levels){auto row=base;row.values[parameter.name]=level;next.push_back(std::move(row));}
        cases=std::move(next);
    }
    for(std::size_t i=0;i<cases.size();++i)cases[i].case_id=case_id(i);
    return cases;
}

std::vector<CampaignCase> generate_latin_hypercube_campaign(std::span<const CampaignParameter> parameters,std::size_t samples,CampaignRandom rng){
    if(parameters.empty()||samples==0U)throw std::invalid_argument("latin hypercube campaign requires parameters and samples");
    for(const auto& p:parameters)validate_parameter(p);
    std::vector<CampaignCase> cases(samples);for(std::size_t i=0;i<samples;++i)cases[i].case_id=case_id(i);
    for(const auto& parameter:parameters){
        if(parameter.kind==CampaignParameterKind::continuous){
            std::vector<double> strata(samples);for(std::size_t i=0;i<samples;++i)strata[i]=(static_cast<double>(i)+rng.uniform_open01())/static_cast<double>(samples);
            for(std::size_t i=samples;i>1U;--i){const std::size_t j=static_cast<std::size_t>(rng.uniform_open01()*static_cast<double>(i));std::swap(strata[i-1U],strata[j]);}
            for(std::size_t i=0;i<samples;++i){const double value=parameter.minimum+(parameter.maximum-parameter.minimum)*strata[i];cases[i].values[parameter.name]=number_to_string(value);} 
        } else {
            for(std::size_t i=0;i<samples;++i)cases[i].values[parameter.name]=parameter.levels[i%parameter.levels.size()];
        }
    }
    return cases;
}

std::string render_campaign_template(std::string_view templ,const CampaignCase& row,bool strict){
    auto replace_placeholders=[&](std::string line){
        std::size_t pos=0U;
        while((pos=line.find('{',pos))!=std::string::npos){const auto end=line.find('}',pos+1U);if(end==std::string::npos)break;const auto key=line.substr(pos+1U,end-pos-1U);auto it=row.values.find(key);if(it==row.values.end()){if(strict)throw std::invalid_argument("missing campaign template placeholder: "+key);pos=end+1U;continue;}line.replace(pos,end-pos+1U,it->second);pos+=it->second.size();}
        return line;
    };
    std::stringstream in{std::string(templ)};std::ostringstream out;std::string line;std::vector<bool> include_stack;
    auto currently_included=[&](){return std::all_of(include_stack.begin(),include_stack.end(),[](bool v){return v;});};
    while(std::getline(in,line)){
        const auto if_pos=line.find("<!-- IF ");
        const auto endif_pos=line.find("<!-- ENDIF -->");
        if(if_pos!=std::string::npos){
            const auto close=line.find("-->",if_pos);if(close==std::string::npos)throw std::invalid_argument("unterminated campaign IF block");std::string expr=line.substr(if_pos+8U,close-(if_pos+8U));
            bool include=false;auto neq=expr.find("!=");auto eq=expr.find('=');
            if(neq!=std::string::npos){const auto key=expr.substr(0,neq);const auto expected=expr.substr(neq+2U);auto it=row.values.find(key);include=it!=row.values.end()&&it->second!=expected;}
            else if(eq!=std::string::npos){const auto key=expr.substr(0,eq);const auto expected=expr.substr(eq+1U);auto it=row.values.find(key);include=it!=row.values.end()&&it->second==expected;}
            else {auto it=row.values.find(expr);include=it!=row.values.end()&&truthy(it->second);}include_stack.push_back(include);continue;
        }
        if(endif_pos!=std::string::npos){if(include_stack.empty())throw std::invalid_argument("campaign ENDIF without IF");include_stack.pop_back();continue;}
        if(currently_included())out<<replace_placeholders(line)<<'\n';
    }
    if(!include_stack.empty())throw std::invalid_argument("unclosed campaign IF block");return out.str();
}

void validate_solver_adapter_descriptor(const SolverAdapterDescriptor& adapter){
    if(adapter.name.empty()||adapter.native_binary.empty()||adapter.container_binary.empty())throw std::invalid_argument("solver adapter descriptor requires identity and binaries");
    for(const auto& action:adapter.control_actions)if(action.empty())throw std::invalid_argument("solver control action must not be empty");
}

bool solver_has_capability(const SolverAdapterDescriptor& adapter,SolverCapability capability){return std::find(adapter.capabilities.begin(),adapter.capabilities.end(),capability)!=adapter.capabilities.end();}

void update_case_status(std::vector<CampaignRegistryEntry>& registry,CampaignRegistryEntry entry){
    if(entry.case_id.empty())throw std::invalid_argument("campaign registry entry needs case_id");
    auto it=std::find_if(registry.begin(),registry.end(),[&](const auto& item){return item.case_id==entry.case_id;});
    if(it==registry.end())registry.push_back(std::move(entry));else *it=std::move(entry);
}

CampaignSummary summarize_campaign_registry(std::span<const CampaignRegistryEntry> registry){
    CampaignSummary out;bool have_best=false;
    for(const auto& item:registry){
        switch(item.status){case CampaignCaseStatus::pending:++out.pending;break;case CampaignCaseStatus::running:++out.running;break;case CampaignCaseStatus::done:++out.done;break;case CampaignCaseStatus::failed:++out.failed;break;case CampaignCaseStatus::stopped:++out.stopped;break;}
        if(item.status==CampaignCaseStatus::done&&(!have_best||item.objective<out.best_objective)){have_best=true;out.best_objective=item.objective;out.best_case_id=item.case_id;}
    }
    return out;
}

std::vector<double> finite_difference_gradient(const std::function<double(std::span<const double>)>& objective,std::span<const double> parameters,const FiniteDifferenceGradientConfig& config){
    if(!objective||parameters.empty()||!(config.absolute_step>0.0)||config.relative_step<0.0)throw std::invalid_argument("invalid finite difference gradient request");
    std::vector<double> x(parameters.begin(),parameters.end()),gradient(x.size());
    const double base=objective(x);
    for(std::size_t i=0;i<x.size();++i){const double h=std::max(config.absolute_step,config.relative_step*std::abs(x[i]));auto xp=x;xp[i]+=h;if(config.central){auto xm=x;xm[i]-=h;gradient[i]=(objective(xp)-objective(xm))/(2.0*h);}else gradient[i]=(objective(xp)-base)/h;}
    return gradient;
}

std::vector<double> gradient_descent_update(std::span<const double> parameters,std::span<const double> gradient,double step_size){
    if(parameters.size()!=gradient.size()||parameters.empty()||!(step_size>=0.0))throw std::invalid_argument("invalid gradient descent update");std::vector<double> out(parameters.begin(),parameters.end());for(std::size_t i=0;i<out.size();++i)out[i]-=step_size*gradient[i];return out;
}

} // namespace cfd::workflow

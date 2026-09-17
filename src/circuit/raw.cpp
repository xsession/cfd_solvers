#include "cfd/circuit/raw.hpp"

#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cfd::circuit {
namespace {
std::string trim(std::string text){
    const auto first=text.find_first_not_of(" \t\r\n");
    if(first==std::string::npos)return {};
    const auto last=text.find_last_not_of(" \t\r\n");
    return text.substr(first,last-first+1U);
}
std::string after_colon(const std::string& line){
    const auto colon=line.find(':');
    if(colon==std::string::npos)throw std::invalid_argument("malformed SPICE RAW header");
    return trim(line.substr(colon+1U));
}
Complex parse_complex(std::string text){
    text=trim(std::move(text));
    if(text.empty())throw std::invalid_argument("empty SPICE RAW value");
    if(text.front()!='(')return {std::stod(text),0.0};
    if(text.back()!=')')throw std::invalid_argument("malformed SPICE RAW complex value");
    const auto comma=text.find(',');
    if(comma==std::string::npos)throw std::invalid_argument("malformed SPICE RAW complex pair");
    return {std::stod(text.substr(1U,comma-1U)),std::stod(text.substr(comma+1U,text.size()-comma-2U))};
}
}

void write_spice_raw(std::ostream& output,const SpiceRawDataset& dataset){
    if(dataset.scale.empty())throw std::invalid_argument("SPICE RAW dataset requires points");
    if(dataset.scale_name.empty()||dataset.scale_type.empty())throw std::invalid_argument("SPICE RAW dataset requires scale metadata");
    for(const auto& trace:dataset.traces){
        if(trace.name.empty()||trace.values.size()!=dataset.scale.size())throw std::invalid_argument("SPICE RAW trace size/name mismatch");
    }
    output<<"Title: "<<dataset.title<<'\n'
          <<"Plotname: "<<dataset.plot_name<<'\n'
          <<"Flags: "<<(dataset.complex_values?"complex":"real")<<'\n'
          <<"No. Variables: "<<dataset.traces.size()+1U<<'\n'
          <<"No. Points: "<<dataset.scale.size()<<'\n'
          <<"Variables:\n"
          <<"\t0\t"<<dataset.scale_name<<'\t'<<dataset.scale_type<<'\n';
    for(std::size_t i=0;i<dataset.traces.size();++i)
        output<<'\t'<<i+1U<<'\t'<<dataset.traces[i].name<<'\t'<<dataset.traces[i].type<<'\n';
    output<<"Values:\n"<<std::setprecision(17);
    for(std::size_t point=0;point<dataset.scale.size();++point){
        if(dataset.complex_values)output<<point<<"\t("<<dataset.scale[point]<<",0)\n";
        else output<<point<<'\t'<<dataset.scale[point]<<'\n';
        for(const auto& trace:dataset.traces){
            const auto value=trace.values[point];
            if(dataset.complex_values)output<<"\t("<<value.real()<<','<<value.imag()<<")\n";
            else output<<'\t'<<value.real()<<'\n';
        }
    }
}

SpiceRawDataset read_spice_raw(std::istream& input){
    SpiceRawDataset dataset;
    std::string line;
    std::size_t variables=0U,points=0U;
    bool in_variables=false,found_values=false;
    while(std::getline(input,line)){
        const std::string clean=trim(line);
        if(clean.empty())continue;
        if(clean.rfind("Title:",0)==0){dataset.title=after_colon(clean);continue;}
        if(clean.rfind("Plotname:",0)==0){dataset.plot_name=after_colon(clean);continue;}
        if(clean.rfind("Flags:",0)==0){dataset.complex_values=after_colon(clean).find("complex")!=std::string::npos;continue;}
        if(clean.rfind("No. Variables:",0)==0){variables=static_cast<std::size_t>(std::stoull(after_colon(clean)));continue;}
        if(clean.rfind("No. Points:",0)==0){points=static_cast<std::size_t>(std::stoull(after_colon(clean)));continue;}
        if(clean=="Variables:"){in_variables=true;continue;}
        if(clean=="Values:"){found_values=true;in_variables=false;break;}
        if(in_variables){
            std::istringstream row(clean);std::size_t index=0U;std::string name,type;
            if(!(row>>index>>name>>type))throw std::invalid_argument("malformed SPICE RAW variable row");
            if(index==0U){dataset.scale_name=name;dataset.scale_type=type;}
            else dataset.traces.push_back({name,type,{}});
        }
    }
    if(!found_values||variables==0U||variables!=dataset.traces.size()+1U||points==0U)
        throw std::invalid_argument("incomplete SPICE RAW header");
    dataset.scale.reserve(points);
    for(auto& trace:dataset.traces)trace.values.reserve(points);
    for(std::size_t point=0;point<points;++point){
        do{if(!std::getline(input,line))throw std::invalid_argument("truncated SPICE RAW values");line=trim(line);}while(line.empty());
        std::istringstream first(line);std::size_t index=0U;std::string value_text;
        if(!(first>>index>>value_text)||index!=point)throw std::invalid_argument("SPICE RAW point index mismatch");
        dataset.scale.push_back(parse_complex(value_text).real());
        for(auto& trace:dataset.traces){
            do{if(!std::getline(input,line))throw std::invalid_argument("truncated SPICE RAW trace values");line=trim(line);}while(line.empty());
            trace.values.push_back(parse_complex(line));
        }
    }
    return dataset;
}

SpiceRawDataset raw_from_transient(std::span<const TransientPoint> points,
                                   std::span<const std::string> node_names,std::string title){
    if(points.empty())throw std::invalid_argument("cannot export empty transient RAW dataset");
    if(points.front().node_voltage.size()!=node_names.size())throw std::invalid_argument("transient RAW node-name mismatch");
    SpiceRawDataset out;out.title=std::move(title);out.plot_name="Transient Analysis";out.scale_name="time";out.scale_type="time";
    out.scale.reserve(points.size());
    for(std::size_t node=1U;node<node_names.size();++node)out.traces.push_back({"v("+node_names[node]+")","voltage",{}});
    for(const auto& point:points){
        if(point.node_voltage.size()!=node_names.size())throw std::invalid_argument("transient RAW topology changed");
        out.scale.push_back(point.time);
        for(std::size_t node=1U;node<node_names.size();++node)out.traces[node-1U].values.push_back({point.node_voltage[node],0.0});
    }
    return out;
}

SpiceRawDataset raw_from_ac(std::span<const AcSolution> points,
                            std::span<const std::string> node_names,std::string title){
    if(points.empty())throw std::invalid_argument("cannot export empty AC RAW dataset");
    if(points.front().node_voltage.size()!=node_names.size())throw std::invalid_argument("AC RAW node-name mismatch");
    SpiceRawDataset out;out.title=std::move(title);out.plot_name="AC Analysis";out.scale_name="frequency";out.scale_type="frequency";out.complex_values=true;
    out.scale.reserve(points.size());
    for(std::size_t node=1U;node<node_names.size();++node)out.traces.push_back({"v("+node_names[node]+")","voltage",{}});
    for(const auto& point:points){
        if(point.node_voltage.size()!=node_names.size())throw std::invalid_argument("AC RAW topology changed");
        out.scale.push_back(point.frequency_hz);
        for(std::size_t node=1U;node<node_names.size();++node)out.traces[node-1U].values.push_back(point.node_voltage[node]);
    }
    return out;
}

} // namespace cfd::circuit

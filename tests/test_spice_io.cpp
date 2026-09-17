#include "cfd/circuit/spice.hpp"
#include "cfd/circuit/raw.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
bool near(double a,double b,double tolerance){return std::abs(a-b)<=tolerance;}
}

int main(){
    try {
        const char* scoped = R"SPICE(
.subckt SWCELL A B CP CN PARAMS: RONV=10
.model SLOCAL SW(RON={RONV} ROFF=1G VT=0.5)
S1 A B CP CN SLOCAL
.ends SWCELL
V1 a 0 1
VC1 c1 0 1
VC0 c0 0 0
X1 a o1 c1 c0 SWCELL RONV=10
R1 o1 0 90
X2 a o2 c1 c0 SWCELL RONV=90
R2 o2 0 10
)SPICE";
        auto scoped_circuit=cfd::circuit::Circuit::parse_spice(scoped);
        const auto scoped_op=scoped_circuit.dc_operating_point();
        require(scoped_op.converged,"scoped-model operating point");
        require(near(scoped_circuit.voltage(scoped_op,"o1"),0.9,1e-9),"first local model scope");
        require(near(scoped_circuit.voltage(scoped_op,"o2"),0.1,1e-9),"second local model scope");

        const auto root=std::filesystem::temp_directory_path()/"cfd_spice_io_test";
        std::filesystem::remove_all(root);
        std::filesystem::create_directories(root);
        {
            std::ofstream params(root/"params.inc");
            params<<".param BASE=1k\n";
        }
        {
            std::ofstream library(root/"parts.lib");
            library<<".lib FAST\n"
                   <<".subckt LOAD A B PARAMS: RLOAD=2k\n"
                   <<"RLOAD A B {RLOAD}\n"
                   <<".ends LOAD\n"
                   <<".endl FAST\n"
                   <<".lib SLOW\n"
                   <<".subckt LOAD A B PARAMS: RLOAD=20k\n"
                   <<"RLOAD A B {RLOAD}\n"
                   <<".ends LOAD\n"
                   <<".endl SLOW\n";
        }
        const auto main_file=root/"main.cir";
        {
            std::ofstream main(main_file);
            main<<".include \"params.inc\"\n"
                <<".lib \"parts.lib\" FAST\n"
                <<".func TWICE(X) {2*X}\n"
                <<"V1 in 0 9\n"
                <<"RUP in out {TWICE(BASE)}\n"
                <<"RDN out 0 BASE\n"
                <<"V2 aux 0 3\n"
                <<"XLOAD aux 0 LOAD RLOAD={3*BASE}\n";
        }
        auto file_circuit=cfd::circuit::Circuit::parse_spice_file(main_file.string());
        const auto file_op=file_circuit.dc_operating_point();
        require(file_op.converged,"file SPICE operating point");
        require(near(file_circuit.voltage(file_op,"out"),3.0,1e-8),".FUNC/.INCLUDE expression result");
        require(near(file_circuit.voltage(file_op,"aux"),3.0,1e-9),".LIB selected subcircuit");

        const char* initial_net=R"SPICE(
R1 out 0 1k
C1 out 0 1u
.IC V(out)=1
.NODESET V(out)=0.25
)SPICE";
        auto initial=cfd::circuit::Circuit::parse_spice(initial_net);
        const auto out_node=initial.find_node("out");
        require(near(initial.initial_voltages().at(out_node),1.0,1e-15),".IC parser stores node initial condition");
        require(near(initial.nodeset_voltages().at(out_node),0.25,1e-15),".NODESET parser stores Newton hint");
        const auto discharge=initial.transient(1e-4,3);
        require(near(discharge.front().node_voltage[out_node],1.0,1e-15)&&discharge.back().node_voltage[out_node]<1.0,
                ".IC initializes transient state");

        const auto raw_dataset=cfd::circuit::raw_from_transient(discharge,initial.node_names());
        std::stringstream raw_text;
        cfd::circuit::write_spice_raw(raw_text,raw_dataset);
        const auto raw_roundtrip=cfd::circuit::read_spice_raw(raw_text);
        require(raw_roundtrip.scale.size()==discharge.size()&&raw_roundtrip.traces.size()==1U,
                "SPICE RAW transient topology roundtrip");
        require(near(raw_roundtrip.traces[0].values.front().real(),1.0,1e-15),"SPICE RAW transient value roundtrip");

        cfd::circuit::Circuit ac_export;const auto acin=ac_export.node("in"),acout=ac_export.node("out");
        ac_export.add_voltage_source("V",acin,0,0.0,{1.0,0.0});ac_export.add_resistor("R",acin,acout,1000.0);ac_export.add_capacitor("C",acout,0,1e-6);
        const auto ac_points=ac_export.ac_log_sweep(10.0,1000.0,5);
        const auto ac_raw=cfd::circuit::raw_from_ac(ac_points,ac_export.node_names());
        std::stringstream ac_text;cfd::circuit::write_spice_raw(ac_text,ac_raw);
        const auto ac_roundtrip=cfd::circuit::read_spice_raw(ac_text);
        require(ac_roundtrip.complex_values&&ac_roundtrip.scale.size()==5U&&ac_roundtrip.traces.size()==2U,
                "SPICE RAW complex AC roundtrip");
        require(std::abs(ac_roundtrip.traces[1].values.back()-ac_points.back().node_voltage[acout])<1e-14,
                "SPICE RAW AC complex value roundtrip");

        const char* behavioral_net=R"SPICE(
.param GAIN=2
V1 in 0 DC 1 AC 1
BGAIN out 0 V={GAIN*V(in)}
RLOAD out 0 1k
BQUAD sense 0 I={1m*V(in)*V(in)}
RSENSE sense 0 1k
BTIME ramp 0 V={1000*TIME}
)SPICE";
        auto behavioral=cfd::circuit::Circuit::parse_spice(behavioral_net);
        const auto behavioral_op=behavioral.dc_operating_point_homotopy();
        require(behavioral_op.converged,"behavioral-source DC operating point");
        require(near(behavioral.voltage(behavioral_op,"out"),2.0,2e-8),"behavioral voltage source DC expression");
        require(near(behavioral.voltage(behavioral_op,"sense"),-1.0,2e-8),"behavioral current source DC expression");
        const auto behavioral_ac=behavioral.ac(1000.0,&behavioral_op);
        require(std::abs(behavioral.voltage(behavioral_ac,"out")-cfd::circuit::Complex{2.0,0.0})<2e-7,
                "behavioral voltage source small-signal linearization");
        require(std::abs(behavioral.voltage(behavioral_ac,"sense")-cfd::circuit::Complex{-2.0,0.0})<2e-6,
                "nonlinear behavioral current source AC Jacobian");
        const auto behavioral_tran=behavioral.transient(1e-3,1U);
        require(near(behavioral_tran.back().node_voltage[behavioral.find_node("ramp")],1.0,2e-8),
                "behavioral TIME expression in transient analysis");

        const char* dependent_nonlinear=R"SPICE(
VIN in 0 DC 2 AC 1
EPOLY poly 0 POLY(1) in 0 1 2 3
ETABLE tab 0 TABLE {V(in)} = (0,0) (1,2) (3,4)
RP poly 0 1k
RT tab 0 1k
)SPICE";
        auto dependent=cfd::circuit::Circuit::parse_spice(dependent_nonlinear);
        const auto dependent_op=dependent.dc_operating_point_homotopy();
        require(dependent_op.converged&&near(dependent.voltage(dependent_op,"poly"),17.0,2e-7),
                "dependent POLY(1) DC syntax");
        require(near(dependent.voltage(dependent_op,"tab"),3.0,2e-7),"dependent TABLE interpolation syntax");
        const auto dependent_ac=dependent.ac(1.0e3,&dependent_op);
        require(std::abs(dependent.voltage(dependent_ac,"poly")-cfd::circuit::Complex{14.0,0.0})<2e-5,
                "dependent POLY(1) small-signal derivative");
        require(std::abs(dependent.voltage(dependent_ac,"tab")-cfd::circuit::Complex{1.0,0.0})<2e-5,
                "dependent TABLE small-signal slope");

        const char* behavioral_subckt=R"SPICE(
.subckt GAINCELL IN OUT PARAMS: K=3
B1 OUT 0 V={K*V(IN)}
.ends GAINCELL
VIN source 0 2
XG source result GAINCELL K=4
RLOAD result 0 1k
)SPICE";
        auto behavioral_hierarchy=cfd::circuit::Circuit::parse_spice(behavioral_subckt);
        const auto hierarchy_op=behavioral_hierarchy.dc_operating_point_homotopy();
        require(hierarchy_op.converged&&near(behavioral_hierarchy.voltage(hierarchy_op,"result"),8.0,1e-7),
                "behavioral source subcircuit node/parameter flattening");

        bool include_rejected=false;
        try { (void)cfd::circuit::Circuit::parse_spice(".include external.cir\n"); }
        catch(const std::invalid_argument&){include_rejected=true;}
        require(include_rejected,"string parser must reject file directives");
        std::filesystem::remove_all(root);
        std::cout<<"SPICE IO tests passed\n";
        return 0;
    } catch(const std::exception& e){
        std::cerr<<"SPICE IO failure: "<<e.what()<<'\n';
        return 1;
    }
}

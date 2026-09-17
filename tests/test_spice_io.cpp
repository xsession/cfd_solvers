#include "cfd/circuit/spice.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
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

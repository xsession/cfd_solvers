#include "cfd/solvers/fdtd/vtk_io.hpp"

#include "cfd/solvers/fdtd/maxwell3d.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace cfd::fdtd {

void write_maxwell3d_vtk_ascii(const Maxwell3D& s,const std::filesystem::path& path){
    std::ofstream out(path);
    if(!out) throw std::runtime_error("failed to open VTK output file: "+path.string());
    out << "# vtk DataFile Version 3.0\n"
        << "cfd_solvers Maxwell3D fields\nASCII\nDATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << s.nx() << ' ' << s.ny() << ' ' << s.nz() << "\n"
        << "ORIGIN 0 0 0\n"
        << std::setprecision(17)
        << "SPACING " << s.dx() << ' ' << s.dy() << ' ' << s.dz() << "\n"
        << "POINT_DATA " << s.cell_count() << "\n"
        << "VECTORS E double\n";
    for(std::size_t n=0;n<s.cell_count();++n) out << s.ex()[n] << ' ' << s.ey()[n] << ' ' << s.ez()[n] << '\n';
    out << "VECTORS H double\n";
    for(std::size_t n=0;n<s.cell_count();++n) out << s.hx()[n] << ' ' << s.hy()[n] << ' ' << s.hz()[n] << '\n';
    if(!out) throw std::runtime_error("failed while writing VTK output file: "+path.string());
}

} // namespace cfd::fdtd

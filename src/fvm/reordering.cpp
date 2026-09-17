#include "cfd/fvm/reordering.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
namespace cfd::fvm {namespace {std::uint64_t part1(std::uint32_t x){std::uint64_t v=x&0x1fffffU;v=(v|(v<<32))&0x1f00000000ffffULL;v=(v|(v<<16))&0x1f0000ff0000ffULL;v=(v|(v<<8))&0x100f00f00f00f00fULL;v=(v|(v<<4))&0x10c30c30c30c30c3ULL;v=(v|(v<<2))&0x1249249249249249ULL;return v;}std::uint64_t code(std::uint32_t x,std::uint32_t y,std::uint32_t z){return part1(x)|(part1(y)<<1)|(part1(z)<<2);}}
std::vector<std::size_t> morton_cell_order(const PolyMesh&m){if(m.cell_count()==0)return {};Vec3 lo=m.cells()[0].center,hi=lo;for(const auto&c:m.cells()){lo.x=std::min(lo.x,c.center.x);lo.y=std::min(lo.y,c.center.y);lo.z=std::min(lo.z,c.center.z);hi.x=std::max(hi.x,c.center.x);hi.y=std::max(hi.y,c.center.y);hi.z=std::max(hi.z,c.center.z);}auto quant=[](double v,double a,double b){if(!(b>a))return 0U;return static_cast<std::uint32_t>(std::clamp((v-a)/(b-a),0.0,1.0)*2097151.0);};std::vector<std::size_t> order(m.cell_count());std::iota(order.begin(),order.end(),0);std::sort(order.begin(),order.end(),[&](auto a,auto b){auto A=m.cells()[a].center,B=m.cells()[b].center;return code(quant(A.x,lo.x,hi.x),quant(A.y,lo.y,hi.y),quant(A.z,lo.z,hi.z))<code(quant(B.x,lo.x,hi.x),quant(B.y,lo.y,hi.y),quant(B.z,lo.z,hi.z));});return order;}
} // namespace cfd::fvm

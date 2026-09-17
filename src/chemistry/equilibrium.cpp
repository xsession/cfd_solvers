#include "cfd/chemistry/equilibrium.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd::chemistry {
BinaryPrecipitationResult equilibrate_binary_salt(double a,double b,double solid,double ksp){
    if(a<0||b<0||solid<0||!(ksp>0)||!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(solid)||!std::isfinite(ksp))throw std::invalid_argument("invalid binary-salt equilibrium");
    const double disc=(a-b)*(a-b)+4.0*ksp;const double equilibrium_extent=0.5*(a+b-std::sqrt(disc));
    double extent=equilibrium_extent;if(extent<0.0)extent=std::max(extent,-solid);extent=std::min(extent,std::min(a,b));
    const double af=a-extent,bf=b-extent,sf=solid+extent;return {af,bf,sf,extent,af*bf/ksp};
}
}

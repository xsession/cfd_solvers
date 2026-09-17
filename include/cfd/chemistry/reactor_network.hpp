#pragma once

#include "cfd/chemistry/implicit_reactor.hpp"
#include <cstddef>
#include <vector>

namespace cfd::chemistry {

struct WellStirredReactor {
    double volume{1.0}; // m3
    double temperature{298.15};
    std::vector<double> concentrations;
};
struct ReactorConnection {
    std::size_t first{};
    std::size_t second{};
    double exchange_flow{}; // m3/s, symmetric mixing exchange
};

class WellStirredReactorNetwork {
public:
    explicit WellStirredReactorNetwork(ReactionNetwork chemistry);
    std::size_t add_reactor(WellStirredReactor reactor);
    void connect(ReactorConnection connection);
    void advance(double duration,double maximum_step=1.0e-3,const ReactorConfig& chemistry_controls={});
    [[nodiscard]] const std::vector<WellStirredReactor>& reactors() const noexcept{return reactors_;}
    [[nodiscard]] std::vector<double> total_species_amounts() const;
private:
    ReactionNetwork chemistry_;
    std::vector<WellStirredReactor> reactors_;
    std::vector<ReactorConnection> connections_;
};

} // namespace cfd::chemistry

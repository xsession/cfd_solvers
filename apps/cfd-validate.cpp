#include "cfd/workflow/validation_catalog.hpp"
#include <iostream>
int main(){std::cout<<"name\tfamily\tcommand\tcriterion\n";for(const auto&c:cfd::workflow::validation_catalog())std::cout<<c.name<<'\t'<<c.family<<'\t'<<c.command<<'\t'<<c.criterion<<'\n';}

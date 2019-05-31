#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "geometry.hpp"

using namespace std;
using namespace arma;



int main()
{
    arma_rng::set_seed(time(NULL));

    // create geometry from input
    Geom24 G(cin);
    
    G.print_omega_table_4();


    return 0;
}

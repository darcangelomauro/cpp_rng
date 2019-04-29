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
    
    // randomize H and L
    G.shuffle();

    // cout dirac op
    cx_mat dirac = G.build_dirac();
    // cout << G.build_dirac() << endl;

    // calculate action

    cout.precision(16);
    cout << G.calculate_S_from_dirac() << endl;
    cout << G.calculate_S() << endl;


    return 0;
}

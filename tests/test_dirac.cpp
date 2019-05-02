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
    Geom24 G(2, 2, 10, 2.2431);
    
    // randomize H and L
    G.shuffle();

    // cout dirac op
    cx_mat dirac = G.build_dirac();
    // cout << G.build_dirac() << endl;

    // calculate action

    cout.precision(16);
    //cout << G.calculate_S_from_dirac() << endl;


    clock_t start1 = clock();
    double S_new = G.calculate_S_new();
    clock_t start2 = clock();
    double S_old = G.calculate_S();
    clock_t end = clock();

    cout << "new: " << S_new << "    time: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "old: " << S_old << "    time: " << (double)(end-start2)/(double)CLOCKS_PER_SEC << endl;

    return 0;
}

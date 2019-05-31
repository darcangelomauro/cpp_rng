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
    Geom24 G(2, 2, 40, 2.2431);
    
    // randomize H and L
    G.shuffle();

    // cout dirac op
    cx_mat dirac = G.build_dirac();
    // cout << G.build_dirac() << endl;

    // calculate action

    cout.precision(16);


    clock_t start1 = clock();
    double S_dir = G.calculate_S_from_dirac();
    clock_t start2 = clock();
    double S_new = G.calculate_S();
    clock_t start3 = clock();
    double S_old = G.calculate_S_old();
    clock_t end = clock();

    cout << "dir: " << S_dir << "    time: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "new: " << S_new << "    time: " << (double)(start3-start2)/(double)CLOCKS_PER_SEC << endl;
    cout << "old: " << S_old << "    time: " << (double)(end-start3)/(double)CLOCKS_PER_SEC << endl;

    return 0;
}

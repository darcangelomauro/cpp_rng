#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;


int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, 222222);
    arma_rng::set_seed(111111);

    // create geometry from input
    Geom24 G(0, 3, 20, -2.5);
    
    ofstream out_da_s;
    out_da_s.open("base_da_s.txt");
    ofstream out_da_hl;
    out_da_hl.open("base_da_hl.txt");
    double dt = 0.001;
    dt = G.dual_averaging(100, dt, 100, engine, out_da_s, out_da_hl);
    out_da_s.close();
    out_da_hl.close();

    cout << "Dual averaging result: " << dt << endl;

    ofstream out_s;
    out_s.open("base_s.txt");
    ofstream out_hl;
    out_hl.open("base_hl.txt");
    double ar = G.HMC(100, dt, 100, engine, out_s, out_hl);
    out_s.close();
    out_hl.close();
    
    cout << "Acceptance rate: " << ar << endl;

    return 0;
}
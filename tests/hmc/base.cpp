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

    // create geometry from input
    Geom24 G(2, 0, 10, -2.5);
    
    ofstream out_s;
    out_s.open("data/base_s.txt");
    ofstream out_hl;
    out_hl.open("data/base_hl.txt");
    double dt = 0.001;
    G.HMC(100, dt, 100, true, engine, out_s, out_hl);
    double ar = G.HMC(100, dt, 100, false, engine, out_s, out_hl);
    out_s.close();
    out_hl.close();

    cout << "acceptance rate: " << ar << endl;


    return 0;
}

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
    gsl_rng_set(engine, time(NULL));

    // create geometry from input
    Geom24 G(2, 2, 20, -2.5);
    
    ofstream out_s;
    out_s.open("data/base_s.txt");

    G.delta24_debug(0.01, 1000, engine, out_s);

    out_s.close();



    return 0;
}

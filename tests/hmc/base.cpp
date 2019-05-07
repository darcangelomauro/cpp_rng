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
    Geom24 G(2, 0, 20, -1.5);
    
    ofstream out;
    out.open("base.txt");
    double ar = G.HMC(1000, 0.0005, 100, engine, out);
    out.close();

    cout << "acceptance rate: " << ar << endl;


    return 0;
}

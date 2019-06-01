#include <iostream>
#include <ostream>
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

    //const double L = 1;

    // create geometry from input
    Geom24 G(1, 1, 20, -2.5);
    
    //G.shuffle();
    G.sample_mom(engine);
    
    ofstream out;
    out.open("data/env.txt");
    out.precision(10);
    double S = G.calculate_S();
    double K = G.calculate_K();
    out << S << " "<< K << " " << S+K << endl;
    for(int i=0; i<100; ++i)
    {
        G.leapfrog(100, 0.0001);
        S = G.calculate_S();
        K = G.calculate_K();
        out << S << " "<< K << " " << S+K << endl;
    }
    out.close();


    return 0;
}

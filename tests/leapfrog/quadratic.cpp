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
    arma_rng::set_seed(time(NULL)+10);

    const double L = 0.01;

    // create geometry from input
    Geom24 G(2, 0, 10, -2.2431);
    
    double tau = 0.00001;

    ofstream out;
    out.open("quadratic.txt");
    out.precision(16);

    while(tau >= 0.000001)
    {
        G.shuffle();
        G.sample_mom(engine);
        
        double Si = G.calculate_S();
        double Ki = G.calculate_K();
        G.leapfrog(L/tau, tau);
        double Sf = G.calculate_S();
        double Kf = G.calculate_K();

        out << tau << " " << fabs(Sf+Kf-Si-Ki) << endl;

        tau /= 1.1;
    }

    out.close();

    return 0;
}

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
    
    double T = 1.;
    double tau = 0.001;

    ofstream out;
    out.open("data/quadratic.txt");
    out.precision(16);

    while(tau >= 0.0001)
    {
        Geom24 G(0, 3, 8, -2.2431);
        G.shuffle(engine);
        G.sample_mom(engine);
        
        /*
        vector<double> err;
        for(int i=0; i<1; ++i)
        {
            double Si = G.calculate_S();
            double Ki = G.calculate_K();
            G.leapfrog(1000, tau);
            double Sf = G.calculate_S();
            double Kf = G.calculate_K();
            err.push_back(fabs(Sf+Kf-Si-Ki));
        }

        out << tau << " " << accumulate(err.begin(), err.end(), 0.)/err.size()  << endl;
        */

        double Si = G.calculate_S();
        double Ki = G.calculate_K();
        G.leapfrog(int(T/tau), tau);
        double Sf = G.calculate_S();
        double Kf = G.calculate_K();
        out << log(tau) << " " << log(fabs(Sf+Kf-Si-Ki)) << endl;
        
        tau /= 1.1;
    }

    out.close();

    return 0;
}

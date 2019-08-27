#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

#define N 1000

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));


    // create geometry
    int p = 2;
    int q = 0;
    int dim = 8;
    double g2 = 2.5;
    
    double tau[4] = {0.001, 0.0005, 0.00025, 0.000125};

    ofstream out;
    string name = filename_from_data(p, q, dim, g2, "LOCERR");
    out.open("data/" + name + ".txt");
    out.precision(16);

    for(int i=0; i<4; ++i)
    {
        Geom24 G(p, q, dim, g2);

        vec dH(N);
        
        for(int j=0; j<N; ++j)
        {
            G.shuffle(engine);
            G.sample_mom(engine);
            
            double Si = G.calculate_S();
            double Ki = G.calculate_K();
            G.leapfrog(1, tau[i]);
            double Sf = G.calculate_S();
            double Kf = G.calculate_K();

            dH(j) = fabs(Sf+Kf-Si-Ki);
        }

        out << log(tau[i]) << " " << log(mean(dH)) << endl;
    }

    out.close();

    return 0;
}

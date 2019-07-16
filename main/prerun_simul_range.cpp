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


int main()
{
    // initialize random number generator
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));





    //********* BEGIN PARAMETER INITIALIZATION **********//

    // geometric parameters
    const int p = 1;
    const int q = 3;
    const int dim = 32;

    // hamiltonian parameters
    const int L = 100;
    double dt = 0.005;
    const int iter_therm = 100;
    
    // coupling constant
    const double g2_i = -3.8;
    const double g2_f = -3.3;
    const double g2_step = 0.1;

    //********* END PARAMETER INITIALIZATION **********//





    string prefix = "HAMIL";
    double g2 = g2_i;
    while(g2 < g2_f)
    {
        Geom24 G(p, q, dim, g2);
       
        cout << "Simulating " << G << endl;

        string path = "data/" + filename_from_data(p, q, dim, g2, prefix);


        // PRELIMINARY RUN

        ofstream out_s, out_hl;
        out_s.open(path + "_S_therm.txt");
        out_hl.open(path + "_HL_therm.txt");

        G.HMC(L, dt, iter_therm, true, engine, out_s, out_hl);
        
        out_s.close();
        out_hl.close();
       
        cout << "dt: " << dt << endl;

        g2 += g2_step;
    }


    return 0;
}

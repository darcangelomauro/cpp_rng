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
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));

    // create geometry from input
    const int p = 1;
    const int q = 3;
    const int dim = 32;
    const double g2_i = -3.8;
    const double g2_f = -3.3;
    const double g2_step = 0.1;
    const int iter_therm = 100;
    const int iter = 1000;

    double g2 = g2_i;
    double dt = 0.005;
    while(g2 < g2_f)
    {
        Geom24 G(p, q, dim, g2);
       
        cout << "Simulating " << G << endl;

        string path = "data/" + filename_from_data(p, q, dim, g2);

        // THERMALIZATION

        ofstream out_s, out_hl;
        out_s.open(path + "_S_therm.txt");
        out_hl.open(path + "_HL_therm.txt");

        G.HMC(100, dt, iter_therm, true, engine, out_s, out_hl);
        
        out_s.close();
        out_hl.close();
       
        thermalization_analysis(path);


        // SIMULATION

        out_s.open(path + "_S.txt");
        out_hl.open(path + "_HL.txt");

        clock_t start1 = clock();
        double ar = G.HMC(100, dt, iter, true, engine, out_s, out_hl);
        clock_t end = clock();

        out_s.close();
        out_hl.close();


        cout << "acceptance rate hmc: " << ar << endl;
        cout << "time hmc: " << (double)(end-start1)/(double)CLOCKS_PER_SEC << endl;

        g2 += g2_step;
    }


    return 0;
}

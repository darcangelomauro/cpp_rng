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
        int c = dim*dim*G.get_nHL();
       
        cout << "Simulating " << G << endl;

        // THERMALIZATION

        ofstream out_s, out_hl;
        out_s.open("data/" + filename_from_data(p, q, dim, g2) + "_S_therm.txt");
        out_hl.open("data/" + filename_from_data(p, q, dim, g2) + "_HL_therm.txt");

        G.HMC(100, dt, iter_therm, true, engine, out_s, out_hl);
        
        out_s.close();
        out_hl.close();
       


        // SIMULATION

        out_s.open("data/" + filename_from_data(p, q, dim, g2) + "_S.txt");
        out_hl.open("data/" + filename_from_data(p, q, dim, g2) + "_HL.txt");

        clock_t start1 = clock();
        double ar = G.HMC(100, dt, iter, true, engine, out_s, out_hl);
        clock_t end = clock();

        out_s.close();
        out_hl.close();

        
        
        // thermalization analysis
        
        double* vec2 = new double [iter_therm];
        double* vec4 = new double [iter_therm];
        for(int i=0; i<iter_therm; ++i)
        {
            vec2[i] = 0;
            vec4[i] = 0;
        }
        
        ifstream in_s;
        in_s.open("data/" + filename_from_data(p, q, dim, g2) + "_S_therm.txt");

        for(int i=0; i<iter_therm; ++i)
            in_s >> vec2[i] >> vec4[i];

        in_s.close();
                


        ofstream out_dofs;
        out_dofs.open("data/" + filename_from_data(p, q, dim, g2) + "_dofs.txt");

        // calculate average of 2gTrD2 + 4TrD4 based on the last iter/10 samples
        for(int i=0; i<(iter_therm-(iter_therm/10)); ++i)
        {
            double res = 0;
            for(int j=0; j<iter_therm/10; ++j)
                res += 2*g2*vec2[i+j] + 4*vec4[i+j];

            out_dofs << 10*res/iter_therm << " " << c << endl;
        }
        out_dofs.close();

        delete [] vec2;
        delete [] vec4;

        cout << "acceptance rate hmc: " << ar << endl;
        cout << "time hmc: " << (double)(end-start1)/(double)CLOCKS_PER_SEC << endl;

        g2 += g2_step;
    }


    return 0;
}

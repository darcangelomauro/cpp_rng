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
    const int p = 2;
    const int q = 0;
    const int dim = 32;
    const double g2 = -2.725;
    const int iter_therm = 100;
    const int iter_sim = 1000;

    Geom24 G1(p, q, dim, g2);
    Geom24 G2(p, q, dim, g2);
    
    int c = dim*dim*G1.get_nHL();
    
    double* vec2_h = new double [iter_therm];
    double* vec4_h = new double [iter_therm];
    double* vec2_m = new double [iter_therm];
    double* vec4_m = new double [iter_therm];
    for(int i=0; i<iter_therm; ++i)
    {
        vec2_h[i] = 0;
        vec4_h[i] = 0;
        vec2_m[i] = 0;
        vec4_m[i] = 0;
    }

    double dt = 0.005;
    double scale = 0.128;


    // THERMALIZATION + SIMULATION

    ofstream out_s_h, out_hl_h, out_s_m, out_hl_m;
    out_s_h.open("data/" + filename_from_data(p, q, dim, g2) + "_S_h.txt");
    out_hl_h.open("data/" + filename_from_data(p, q, dim, g2) + "_HL_h.txt");
    out_s_m.open("data/" + filename_from_data(p, q, dim, g2) + "_S_m.txt");
    out_hl_m.open("data/" + filename_from_data(p, q, dim, g2) + "_HL_m.txt");
    
    G1.HMC_analytic_test(100, dt, iter_therm, true, engine, vec2_h, vec4_h);
    G2.MMC_analytic_test(scale, iter_therm, true, engine, vec2_m, vec4_m);
    
    clock_t start1 = clock();
    double ar_h = G1.HMC(100, dt, iter_sim, false, engine, out_s_h, out_hl_h);
    clock_t start2 = clock();
    double ar_m = G2.MMC(scale, iter_sim, false, engine, out_s_m, out_hl_m);
    clock_t end = clock();
    
    out_s_h.close();
    out_hl_h.close();
    out_s_m.close();
    out_hl_m.close();



    // OUTPUT THERMALIZATION DATA

    ofstream out_s;
    out_s.open("data/" + filename_from_data(p, q, dim, g2) + "_therm.txt");

    // calculate average of 2gTrD2 + 4TrD4 based on the last iter/10 samples
    for(int i=0; i<(iter_therm-(iter_therm/10)); ++i)
    {
        double res_h = 0;
        double res_m = 0;
        for(int j=0; j<iter_therm/10; ++j)
        {
            res_h += 2*g2*vec2_h[i+j] + 4*vec4_h[i+j];
            res_m += 2*g2*vec2_m[i+j] + 4*vec4_m[i+j];

        }
        out_s << 10*res_h/iter_therm << " " << 10*res_m/iter_therm << " " << c << endl;
    }
    out_s.close();

    delete [] vec2_h;
    delete [] vec4_h;
    delete [] vec2_m;
    delete [] vec4_m;

    cout << "hmc dt: " << dt << endl;
    cout << "acceptance rate hmc: " << ar_h << endl;
    cout << "time hmc: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "hmm scale: " << scale << endl;
    cout << "acceptance rate mmc: " << ar_m << endl;
    cout << "time mmc: " << (double)(end-start2)/(double)CLOCKS_PER_SEC << endl;


    return 0;
}

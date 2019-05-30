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
    arma_rng::set_seed(time(NULL)+2);

    // create geometry from input
    const int p = 2;
    const int q = 0;
    const int dim = 30;
    const double g2 = -2.5;
    const int iter = 100;

    Geom24 G1(p, q, dim, g2);
    Geom24 G2(p, q, dim, g2);
    int c = dim*dim*G1.get_nHL();
    
    double* vec2_h = new double [iter];
    double* vec4_h = new double [iter];
    double* vec2_m = new double [iter];
    double* vec4_m = new double [iter];
    for(int i=0; i<iter; ++i)
    {
        vec2_h[i] = 0;
        vec4_h[i] = 0;
        vec2_m[i] = 0;
        vec4_m[i] = 0;
    }

    clock_t start1 = clock();
    double ar_h = G1.HMC_analytic_test(100, 0.0014, iter, engine, vec2_h, vec4_h);
    clock_t start2 = clock();
    double ar_m = G2.MMC_analytic_test(0.041, iter, engine, vec2_m, vec4_m);
    clock_t end = clock();
    
    ofstream out_s;
    out_s.open("p2q0dim10gn2_5.txt");

    // calculate average of 2gTrD2 + 4TrD4 based on the last iter/10 samples
    for(int i=0; i<(iter-(iter/10)); ++i)
    {
        double res_h = 0;
        double res_m = 0;
        for(int j=0; j<iter/10; ++j)
        {
            res_h += 2*g2*vec2_h[i+j] + 4*vec4_h[i+j];
            res_m += 2*g2*vec2_m[i+j] + 4*vec4_m[i+j];

        }
        out_s << 10*res_h/iter << " " << 10*res_m/iter << " " << c << endl;
    }
    out_s.close();

    delete [] vec2_h;
    delete [] vec4_h;
    delete [] vec2_m;
    delete [] vec4_m;

    cout << "acceptance rate hmc: " << ar_h << endl;
    cout << "time hmc: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "acceptance rate mmc: " << ar_m << endl;
    cout << "time mmc: " << (double)(end-start2)/(double)CLOCKS_PER_SEC << endl;


    return 0;
}

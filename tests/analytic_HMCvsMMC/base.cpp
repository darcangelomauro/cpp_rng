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
    const int dim = 20;
    const double g2 = -4.0;
    const int iter = 10000;

    Geom24 G1(p, q, dim, g2);
    Geom24 G2(p, q, dim, g2);
    int c = dim*dim*G1.get_nHL();
   

    // dual averaging hmc
    double dt = 0.005;
    G1.HMC(100, dt, 1000, engine, 0.65);
    
    // dual averaging mmc
    double scale = 0.128;
    G2.MMC(scale, 1000, engine, 0.232);

    ofstream out_s_h, out_s_m;
    out_s_h.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_S.txt");
    out_s_m.open("data/" + filename_from_data(p, q, dim, g2, "MMC") + "_S.txt");


    clock_t start1 = clock();
    double ar_h = G1.HMC(100, dt, iter, 1, engine, out_s_h);
    clock_t start2 = clock();
    double ar_m = G2.MMC(scale, iter, 1, engine, out_s_m);
    clock_t end = clock();
    
    out_s_h.close();
    out_s_m.close();




    // thermalization analysis
    
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
    
    ifstream in_s_h, in_s_m;
    in_s_h.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_S.txt");
    in_s_m.open("data/" + filename_from_data(p, q, dim, g2, "MMC") + "_S.txt");

    for(int i=0; i<iter; ++i)
    {
        in_s_h >> vec2_h[i] >> vec4_h[i];
        in_s_m >> vec2_m[i] >> vec4_m[i];
    }

    in_s_h.close();
    in_s_m.close();
            


    ofstream out_s;
    out_s.open("data/" + filename_from_data(p, q, dim, g2, "") + "_dofs.txt");

    // calculate progression in the average of 2gTrD2 + 4TrD4
    for(int i=0; i<iter; i=i+10)
    {
        double res_h = 0;
        double res_m = 0;
        for(int j=0; j<i; ++j)
        {
            res_h += 2*g2*vec2_h[j] + 4*vec4_h[j];
            res_m += 2*g2*vec2_m[j] + 4*vec4_m[j];

        }
        out_s << res_h/i << " " << res_m/i << " " << c << endl;
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

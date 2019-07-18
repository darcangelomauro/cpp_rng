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
    
    double dt = 0.005;
    double scale = 0.128;


    // THERMALIZATION

    ofstream out_s_h, out_hl_h, out_s_m, out_hl_m;
    out_s_h.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_S_therm.txt");
    out_s_m.open("data/" + filename_from_data(p, q, dim, g2, "MMC") + "_S_therm.txt");
    
    G1.HMC(100, dt, iter_therm, engine, out_s_h);
    G2.MMC(scale, iter_therm, engine, out_s_m);

    out_s_h.close();
    out_hl_h.close();
    out_s_m.close();
    out_hl_m.close();

    
    // SIMULATION
    
    out_s_h.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_S.txt");
    out_hl_h.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_HL.txt");
    out_s_m.open("data/" + filename_from_data(p, q, dim, g2, "MMC") + "_S.txt");
    out_hl_m.open("data/" + filename_from_data(p, q, dim, g2, "MMC") + "_HL.txt");
    
    clock_t start1 = clock();
    double ar_h = G1.HMC(100, dt, iter_sim, 1, engine, out_s_h, out_hl_h);
    clock_t start2 = clock();
    double ar_m = G2.MMC(scale, iter_sim, 1, engine, out_s_m, out_hl_m);
    clock_t end = clock();
    
    out_s_h.close();
    out_hl_h.close();
    out_s_m.close();
    out_hl_m.close();


    cout << "hmc dt: " << dt << endl;
    cout << "acceptance rate hmc: " << ar_h << endl;
    cout << "time hmc: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "hmm scale: " << scale << endl;
    cout << "acceptance rate mmc: " << ar_m << endl;
    cout << "time mmc: " << (double)(end-start2)/(double)CLOCKS_PER_SEC << endl;


    return 0;
}

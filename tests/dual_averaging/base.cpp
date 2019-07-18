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
    gsl_rng_set(engine, 222222);

    const int p = 3;
    const int q = 0;
    const int dim = 8;
    const double g2 = -2.5;

    // create geometry from input
    Geom24 G1(p, q, dim, g2);
    
    ofstream out_da_s_h;
    out_da_s_h.open("data/base_da_s_h.txt");
    double dt = 0.05;
    G1.HMC(100, dt, 100, engine, out_da_s_h);
    out_da_s_h.close();

    cout << "Dual averaging result hmc: " << dt << endl;

    ofstream out_s_h;
    out_s_h.open("data/base_s_h.txt");
    ofstream out_hl_h;
    out_hl_h.open("data/base_hl_h.txt");
    double ar = G1.HMC(100, dt, 100, 1, engine, out_s_h, out_hl_h);
    out_s_h.close();
    out_hl_h.close();
    
    cout << "Acceptance rate hmc: " << ar << endl;
    
    // create geometry from input
    Geom24 G2(p, q, dim, g2);
    
    ofstream out_da_s_m;
    out_da_s_m.open("data/base_da_s_m.txt");
    double scale = 0.05;
    G2.MMC(scale, 100, engine, out_da_s_m);
    out_da_s_m.close();

    cout << "Dual averaging result mmc: " << scale << endl;

    ofstream out_s_m;
    out_s_m.open("data/base_s_m.txt");
    ofstream out_hl_m;
    out_hl_m.open("data/base_hl_m.txt");
    ar = G2.MMC(scale, 100, 1, engine, out_s_m, out_hl_m);
    out_s_m.close();
    out_hl_m.close();
    
    cout << "Acceptance rate mmc: " << ar << endl;
    return 0;
}

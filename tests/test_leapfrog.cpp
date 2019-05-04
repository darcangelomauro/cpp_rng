#include <iostream>
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
    gsl_rng* engine;
    engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, 222222);
    arma_rng::set_seed(111111);

    // create geometry from input
    Geom24 G(2, 2, 10, -2.2431);
    
    G.shuffle();
    G.sample_mom(engine);
    
    double Si = G.calculate_S();
    double Ki = G.calculate_K();
    clock_t start1 = clock();
    G.leapfrog(1, 0.01);
    clock_t end = clock();
    double Sf = G.calculate_S();
    double Kf = G.calculate_K();
    cout << "time: " << (double)(end-start1)/(double)CLOCKS_PER_SEC << endl;
    cout << "initial potential|kinetic|energy: " << Si << "|" << Ki << "|" << Si+Ki << endl;
    cout << "final   potential|kinetic|energy: " << Sf << "|" << Kf << "|" << Sf+Kf << endl;


    return 0;
}

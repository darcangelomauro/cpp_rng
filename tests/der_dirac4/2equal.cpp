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

#define N 1

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));


    // create geometry
    int p = 0;
    int q = 3;
    int dim = 8;
    double g2 = 2.5;
    
    ofstream out;
    string name = filename_from_data(p, q, dim, g2, "2EQUAL");
    out.open("data/" + name + ".txt");
    out.precision(16);

    Geom24 G(p, q, dim, g2);

    G.shuffle(engine);
    G.sample_mom(engine);

    cx_mat M1 = G.debug_2equal_explicit(0);
    cx_mat M2 = G.debug_2equal(0);

    out << M1 << endl << endl;
    out << M2 << endl << endl;


    out.close();

    return 0;
}

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

#define N 10
#define tol 1e-8

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));


    // create geometry
    Geom24 G(cin);
    int p = G.get_p();
    int q = G.get_q();
    int dim = G.get_dim();
    double g2 = G.get_g2();
    
    ofstream out;
    string name = filename_from_data(p, q, dim, g2, "DER4");
    out.open("data/" + name + ".txt");
    out.precision(16);

    for(int i=0; i<N; ++i)
    {
        G.shuffle(engine);
        G.sample_mom(engine);
        
        for(int j=0; j<G.get_nHL(); ++j)
        {
            cx_mat M1 = G.der_dirac4(j, true);
            cx_mat M2 = G.der_dirac4_bruteforce(j, true);
            cx_mat M3 = G.der_dirac4(j, false);
            cx_mat M4 = G.der_dirac4_bruteforce(j, false);
            out << approx_equal(M1, M2, "absdiff", tol) << " " << approx_equal(M3, M4, "absdiff", tol) << " ";
            out << approx_equal(M1, M3, "absdiff", tol) << " " << approx_equal(M2, M4, "absdiff", tol) << " ";
            out << approx_equal(M1, M4, "absdiff", tol) << " " << approx_equal(M2, M3, "absdiff", tol) << endl;
        }

    }


    out.close();

    return 0;
}

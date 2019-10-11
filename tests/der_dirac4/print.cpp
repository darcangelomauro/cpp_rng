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


    // create geometry
    cout << "input p, q, dim, g2" << endl;
    Geom24 G(cin);
    
    G.shuffle(engine);
    G.sample_mom(engine);
    
    cout << "input matrix index" << endl;
    int j;
    cin >> j;
    cout.precision(16);
    cout << scientific;
    G.der_dirac4_cout(j);
    G.der_dirac2_cout(j);

    return 0;
}

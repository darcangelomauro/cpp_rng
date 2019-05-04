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
    gsl_rng_set(engine, 918273);

    // create geometry from input
    Geom24 G(2, 2, 40, 2.2431);
    
    // sample momentum
    clock_t start1 = clock();
    for(int i=0; i<1000; i++)
        G.sample_mom(engine);
    clock_t end = clock();
    cout << "time: " << (double)(end-start1)/(double)CLOCKS_PER_SEC << endl;

    for(int i=0; i<G.get_nHL(); ++i)
        cout << "mom[" << i << "] hermitian? " << G.get_mom(i).is_hermitian() << endl;

    return 0;
}

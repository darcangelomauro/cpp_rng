#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "geometry.hpp"

using namespace std;
using namespace arma;


int main()
{
    arma_rng::set_seed(time(NULL));
    
    const double diff = 0.000001;
    
    // create geometry from input
    Geom24 G(2, 2, 60, 2.2431);
    
    // randomize H and L
    G.shuffle();

    cout.precision(16);
    for(int i=0; i<G.get_nHL(); ++i)
    {
        clock_t start1 = clock();
        cx_mat der0 = G.der_dirac4(i, 0);
        clock_t start2 = clock();
        cx_mat der1 = G.der_dirac4(i, 1);
        clock_t end = clock();
        cout << "dD4/dM" << i << " (0): " << der0.is_hermitian() << "    time: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl;
        cout << "dD4/dM" << i << " (1): " << der1.is_hermitian() << "    time: " << (double)(end-start2)/(double)CLOCKS_PER_SEC << endl;
        cout << "equal up to " << diff << "?: " << approx_equal(der0, der1, "absdiff", diff) << endl; 
    }
    for(int i=0; i<G.get_nHL(); ++i)
    {
        clock_t start1 = clock();
        cx_mat der2 = G.der_dirac2(i);
        clock_t end = clock();
        cout << "dD2/dM " << i << ": " << der2.is_hermitian() << "    time: " << (double)(end-start1)/(double)CLOCKS_PER_SEC << endl;
    }

    return 0;
}

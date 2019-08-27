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
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));

    // create geometry from input
    Geom24 G(0, 3, 8, -2.2431);
    
    //G.shuffle();
    G.sample_mom(engine);
    Geom24 G_old = G;
    
    double Si = G.calculate_S();
    double Ki = G.calculate_K();
    G.leapfrog(1000, 0.0001);
    G.reverse_mom(); 
    G.leapfrog(1000, 0.0001);
    G.reverse_mom(); 
    double Sf = G.calculate_S();
    double Kf = G.calculate_K();
    cout.precision(10);
    cout << "initial potential|kinetic|energy: " << Si << "|" << Ki << "|" << Si+Ki << endl;
    cout << "final   potential|kinetic|energy: " << Sf << "|" << Kf << "|" << Sf+Kf << endl;
    for(int i=0; i<G.get_nHL(); ++i)
    {
        cout << "mat[" + to_string(i) + "]: " << approx_equal(G_old.get_mat(i), G.get_mat(i), "absdiff", 0.0001) << endl;
        cout << "mom[" + to_string(i) + "]: " << approx_equal(G_old.get_mom(i), G.get_mom(i), "absdiff", 0.0001) << endl;
    }



    return 0;
}

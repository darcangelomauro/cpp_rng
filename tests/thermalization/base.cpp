#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"

using namespace std;
using namespace arma;


int main(int argc, char** argv)
{
    if(argc < 2)
    {
        cout << "Too few arguments. Needs filename. Quitting" << endl;
        return 1;
    }

    string basename(argv[1]);

    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));

    // create geometry from input
    Geom24 G(1, 3, 20, -2.5);
    G.shuffle(engine);
    
    // action output file
    ofstream out_s;
    out_s.open("data/" + basename + "_simS.txt");
    out_s.precision(15);
    
    // mat output file
    ofstream out_hl;
    out_hl.open("data/" + basename + "_simHL.txt");
    out_hl.precision(15);
    
    // simulation
    double dt = 0.004;
    G.HMC(100, dt, 100, true, engine, out_s, out_hl);
    double ar = G.HMC(100, dt, 100, false, engine, out_s, out_hl);
    
    out_s.close();
    out_hl.close();

    cout << "acceptance rate: " << ar << endl;


    return 0;
}

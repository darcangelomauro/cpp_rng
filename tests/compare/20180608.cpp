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


int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));
    arma_rng::set_seed(time(NULL)+10);

    double g = -1.5;
    const double step = -0.15;

    ofstream out_1;
    out_1.open("20180608_data.txt");

    for(int i=0; i<15; i++)
    {

        // create geometry from input
        Geom24 G(2, 0, 20, g);
        out_1 << i << " " << G.get_p() << " " << G.get_q() << " " << G.get_dim() << " " << G.get_g2() << endl;
        
        // action output file
        ofstream out_s;
        string filename_s = "20180608_simS_" + to_string(i) + ".txt"; 
        out_s.open(filename_s);
        out_s.precision(15);
        
        // mat output file
        ofstream out_hl;
        string filename_hl = "20180608_simHL_" + to_string(i) + ".txt"; 
        out_hl.open(filename_hl);
        out_hl.precision(15);
        
        // simulation
        double ar = G.HMC(100, 0.0001, 100, engine, out_s, out_hl);
        
        out_s.close();
        out_hl.close();

        cout << "acceptance rate: " << ar << endl;

        g += step;
    }

    out_1.close();


    return 0;
}

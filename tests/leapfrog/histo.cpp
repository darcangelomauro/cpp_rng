#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"
#include "statistics.hpp"

using namespace std;
using namespace arma;

#define N 10000
#define nbins 20

int main()
{
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(engine, time(NULL));


    // create geometry
    string prefix;
    cout << "insert p, q, dim, g2" << endl;
    Geom24 G(cin);
    cout << "insert prefix" << endl;
    cin >> prefix;
    int p = G.get_p();
    int q = G.get_q();
    int dim = G.get_dim();
    double g2 = G.get_g2();
    
    double tau = 0.001;

    ofstream out;
    string name = filename_from_data(p, q, dim, g2, prefix);
    out.open("data/" + name + ".txt");
    out.precision(16);

    G.shuffle(engine);
    G.sample_mom(engine);
    
    vec dH(N);
    double Si = G.calculate_S();
    double Ki = G.calculate_K();
    
    for(int j=0; j<N; ++j)
    {
        Geom24 G1 = G;
        
        G1.leapfrog(100, tau);
        double Sf = G.calculate_S();
        double Kf = G.calculate_K();

        dH(j) = Sf+Kf-Si-Ki;
    }
    
    vec dH_sorted = sort(dH);
    uvec count = hist(dH_sorted, nbins);

    double range = dH_sorted(N-1) - dH_sorted(0);
    double step = range / nbins;

    double center = dH_sorted(0) + step/2;
    for(int i=0; i<nbins; ++i)
    {
        out << center << " " << count(i) << endl;
        center += step;
    }

    out.close();

    return 0;
}

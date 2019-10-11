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

#define N 100
#define M 10

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
    
    double time = 0.0001;
    double tau = 0.00001;

    ofstream out;
    string name = filename_from_data(p, q, dim, g2, prefix);
    out.open("data/" + name + ".txt");
    out.precision(16);

    for(int i=0; i<M; ++i)
    {
        vec dH(N);
        
        for(int j=0; j<N; ++j)
        {
            G.shuffle(engine);
            G.sample_mom(engine);
            
            double Si = G.calculate_S();
            double Ki = G.calculate_K();
            G.leapfrog(int(time/tau), tau);
            double Sf = G.calculate_S();
            double Kf = G.calculate_K();

            dH(j) = log(fabs(Sf+Kf-Si-Ki));
        }

        double avg, var;
        jackknife(dH, avg, var, my_mean);
        double err = sqrt(var);

        out << log(tau) << " " << avg << " " << err << endl;
        cout << tau << " " << exp(avg) << endl;
        
        tau /= 2;
    }

    out.close();

    return 0;
}

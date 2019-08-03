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

    // create geometry from input
    const int p = 2;
    const int q = 0;
    const int dim = 20;
    const double g2 = -4.0;
    const int iter = 1000;
    const int iter_t = 1000;
    const int M = 10;


    Geom24 G(p, q, dim, g2);
    Geom24 G_fn(p, q, dim, g2);
    Geom24 G_fs(p, q, dim, g2);
    Geom24 G_rn(p, q, dim, g2);
    Geom24 G_rs(p, q, dim, g2);
    int c = dim*dim*G.get_nHL();
   

    // dual averaging + thermalization
    double dt_nosplit = 0.005;
    G_fn.HMC_fix_nosplit(100, dt_nosplit, iter_t, engine, 0.65);
    G_fn.HMC_fix_nosplit(100, dt_nosplit, iter_t, engine);
    
    double dt_split = 0.005;
    G_fs.HMC_fix_split(100, dt_split, M, iter_t, engine, 0.65);
    G_fs.HMC_fix_split(100, dt_split, M, iter_t, engine);
    
    double dt_min_nosplit = dt_nosplit;
    G_rn.HMC_fix_nosplit(100, dt_min_nosplit, iter_t, engine, 0.9);
    
    double dt_max_nosplit = dt_nosplit;
    G_rn.HMC_fix_nosplit(100, dt_max_nosplit, iter_t, engine, 0.4);
    
    double dt_min_split = dt_split;
    G_rs.HMC_fix_split(100, dt_min_split, M, iter_t, engine, 0.9);
    
    double dt_max_split = dt_split;
    G_rs.HMC_fix_split(100, dt_max_split, M, iter_t, engine, 0.4);

    
    ofstream out_s_fn;
    ofstream out_s_fs;
    ofstream out_s_rn;
    ofstream out_s_rs;
    out_s_fn.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_fix_nosplit_S.txt");
    out_s_fs.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_fix_split_S.txt");
    out_s_rn.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_rand_nosplit_S.txt");
    out_s_rs.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_rand_split_S.txt");


    clock_t start1 = clock();
    double ar_fn = G_fn.HMC_fix_nosplit(100, dt_nosplit, iter, 1, engine, out_s_fn);
    clock_t start2 = clock();
    double ar_fs = G_fs.HMC_fix_split(100, dt_split, M, iter, 1, engine, out_s_fs);
    clock_t start3 = clock();
    double ar_rn = G_rn.HMC_rand_nosplit(100, 100, dt_min_nosplit, dt_max_nosplit, iter, 1, engine, out_s_rn);
    clock_t start4 = clock();
    double ar_rs = G_rs.HMC_rand_split(100, 100, dt_min_split, dt_max_split, M, iter, 1, engine, out_s_rs);
    clock_t start5 = clock();
    
    out_s_fn.close();
    out_s_fs.close();
    out_s_rn.close();
    out_s_rs.close();

    cout << "acceptance rate hmc_fn: " << ar_fn << endl;
    cout << "integration step hmc_fn: " << dt_nosplit << endl;
    cout << "time hmc_fn: " << (double)(start2-start1)/(double)CLOCKS_PER_SEC << endl << endl;
    
    cout << "acceptance rate hmc_fs: " << ar_fs << endl;
    cout << "integration step hmc_fs: " << dt_split << endl;
    cout << "time hmc_fs: " << (double)(start3-start2)/(double)CLOCKS_PER_SEC << endl << endl;
    
    cout << "acceptance rate hmc_rn: " << ar_rn << endl;
    cout << "integration step hmc_rn: " << dt_min_nosplit << " " << dt_max_nosplit << endl;
    cout << "time hmc_rn: " << (double)(start4-start3)/(double)CLOCKS_PER_SEC << endl << endl;
    
    cout << "acceptance rate hmc_rs: " << ar_rs << endl;
    cout << "integration step hmc_rs: " << dt_min_split << " " << dt_max_split << endl;
    cout << "time hmc_rs: " << (double)(start5-start4)/(double)CLOCKS_PER_SEC << endl << endl;



    // dofs analysis
    string names[4] = {"fix_nosplit", "fix_split", "rand_nosplit", "rand_split"};

    for(int idx=0; idx<4; ++idx)
    {
        cout << "DOFS processing: " << names[idx] << endl;

        double* vec2 = new double [iter];
        double* vec4 = new double [iter];
        for(int i=0; i<iter; ++i)
        {
            vec2[i] = 0;
            vec4[i] = 0;
        }
        
        ifstream in_s;
        in_s.open("data/" + filename_from_data(p, q, dim, g2, "HMC") + "_" + names[idx] + "_S.txt");

        for(int i=0; i<iter; ++i)
            in_s >> vec2[i] >> vec4[i];

        in_s.close();
                
        
        ofstream out_s;
        out_s.open("data/" + filename_from_data(p, q, dim, g2, "") + "_" + names[idx] + "_dofs.txt");

        // calculate progression in the average of 2gTrD2 + 4TrD4
        for(int i=1; i<iter; i=i+10)
        {
            double res = 0;
            for(int j=0; j<i; ++j)
            {
                res += 2*g2*vec2[j] + 4*vec4[j];
            }
            out_s << res/i << " " << c << endl;
        }
        out_s.close();

        delete [] vec2;
        delete [] vec4;
    }



    return 0;
}

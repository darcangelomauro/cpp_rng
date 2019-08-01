#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "geometry.hpp"
#include "clifford.hpp"

using namespace std;
using namespace arma;


double Geom24::HMC_RP_core(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, gsl_rng* engine, double* en_i, double* en_f)
{
    // acceptance probability (return value)
    double e;
    
    // resample momentum
    sample_mom(engine);

    /*
    // sample stepsize uniformly from [dt_min, dt_max)
    double dt = dt_min + (dt_max-dt_min)*gsl_rng_uniform(engine);
    */

    // choose randomly dt_min or dt_max
    double dt;
    double r = gsl_rng_uniform(engine);
    if(r>0.5)
        dt = dt_max;
    else
        dt = dt_min;

    // choose uniformly from [Nt_min, Nt_max)
    double Nt = Nt_min + (Nt_max-Nt_min)*gsl_rng_uniform(engine);


    // store previous configuration
    cx_mat* mat_bk = new cx_mat [nHL];
    for(int j=0; j<nHL; j++)
        mat_bk[j] = mat[j];

    // calculate initial hamiltonian
    en_i[2] = calculate_K();
    en_i[3] = g2*en_i[0]+en_i[1]+en_i[2];

    // leapfrog
    leapfrog(Nt, dt);

    // calculate final hamiltonian
    en_f[0] = dirac2();
    en_f[1] = dirac4();
    en_f[2] = calculate_K();
    en_f[3] = g2*en_f[0]+en_f[1]+en_f[2];


    // metropolis test
    
    // sometimes leapfrog diverges and Hf becomes nan.
    // so first of all address this case
    if(std::isnan(en_f[3]))
    {
        e = 0;
        // restore old configuration
        for(int j=0; j<nHL; ++j)
        {
            mat[j] = mat_bk[j];
            en_f[0] = en_i[0];
            en_f[1] = en_i[1];
            en_f[2] = en_i[2];
            en_f[3] = en_i[3];
        }
    }
    // now do the standard metropolis test
    else if(en_f[3] > en_i[3])
    {
        double r = gsl_rng_uniform(engine);
        e = exp(en_i[3]-en_f[3]);

        if(r > e)
        {
            // restore old configuration
            for(int j=0; j<nHL; ++j)
            {
                mat[j] = mat_bk[j];
                en_f[0] = en_i[0];
                en_f[1] = en_i[1];
                en_f[2] = en_i[2];
                en_f[3] = en_i[3];
            }
        }
    }
    else
        e = 1;

    delete [] mat_bk;

    return e;
}

// HMC routine that doesn't performs dual averaging and outputs S2, S4, H, L
double Geom24::HMC_RP(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, const int& gap, gsl_rng* engine, ostream& out_s, ostream& out_hl)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_RP_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
        
        // print once every "gap" iterations
        if( !(i%gap) )
        {
            // print S2 and S4
            out_s << en_f[0] << " " << en_f[1] << endl;
            
            // print mat
            for(int j=0; j<nHL; ++j)
            {
                for(int k=0; k<dim; ++k)
                {
                    for(int l=0; l<dim; ++l)
                        out_hl << mat[j](k,l).real() << " " << mat[j](k,l).imag() << " ";
                }
                out_hl << endl;
            }
        }
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}

// HMC routine that doesn't performs dual averaging and outputs S2, S4
double Geom24::HMC_RP(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, const int& gap, gsl_rng* engine, ostream& out_s)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_RP_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
        
        // print once every "gap" iterations
        if( !(i%gap) )
        {
            // print S2 and S4
            out_s << en_f[0] << " " << en_f[1] << endl;
        }
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}

// HMC routine that doesn't performs dual averaging and doesn't output
double Geom24::HMC_RP(const int& Nt_min, const int& Nt_max, const double& dt_min, const double& dt_max, const int& iter, gsl_rng* engine)
{
    // initial (_i) and final (_f) potential2, potential4, kinetic, hamiltonian 
    double* en_i = new double [4];
    double* en_f = new double [4];

    // return statistic
    double Stat = 0;
    
    // iter repetitions of leapfrog
    for(int i=0; i<iter; ++i)
    {
        // if it's not the first interation set potential to
        // previous final value, otherwise compute it
        if(i)
        {
            en_i[0] = en_f[0];
            en_i[1] = en_f[1];
        }
        else
        {
            en_i[0] = dirac2();
            en_i[1] = dirac4();
        }

        
        // core part of HMC
        Stat += HMC_RP_core(Nt_min, Nt_max, dt_min, dt_max, engine, en_i, en_f);
    }

    delete [] en_i;
    delete [] en_f;

    return (Stat/iter);
}

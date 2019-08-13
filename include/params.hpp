#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <iostream>
#include <string>

// Struct that packs the simulation
// parameters for hmc all together.
// 
// Explanation of integration parameters:
//
// L    =   # of integration steps
// 
// dL   =   randomize L in interval [L-dL, L+dL]
//
// AR   =   target acceptance rate
// 
// dAR  =   randomize dt in interval [dt_min, dt_max]
//          such that dt_min gives AR+dAR and
//          dt_max gives AR-dAR
//
// M    =   split hamiltonian in M iterations of
//          S2 for each iteration of S4
//
// Explanation of other parameters:
//
// iter_therm   =   # of thermalization iterations
// 
// iter_simul   =   # of data-collecting iterations
//
// gap          =   # of iterations to be skipped between
//                  two measurements
//
// adj          =   # of iterations to be skipped between
//                  two applications of hermitization + tracelessization
//
// Explanation of HMC mode:
//
// fix_nosplit  =   don't randomize dt or L, don't split hamiltonian
//
// fix_split    =   don't randomize dt or L, split hamiltonian
//
// rand_nosplit =   randomize dt and L, don't split hamiltonian
//
// rand_split   =   randomize dt and L, split hamiltonian
//
struct Simul_params
{
    // Geometric parameters
    int p;
    int q;
    int dim;

    // Integration parameters
    int L;
    int dL;
    double AR;
    double dAR;
    int M;

    // Other parameters
    int iter_therm;
    int iter_simul;
    int gap;
    int adj;

    double g2_i;
    double g2_f;
    double g2_step;

    // HMC mode
    std::string mode;

    // Validity string
    std::string valid;
    std::string control = "p:q:dim:L:dL:AR:dAR:M:iter_therm:iter_simul:gap:adj:g2_i:g2_f:g2_step:mode:";
};

// Function to read simulation parameters from stream
bool read_init_stream(std::istream&, struct Simul_params&);

// Function to check whether parameters have been entered in the
// correct order
bool params_validity(struct Simul_params&);


#endif


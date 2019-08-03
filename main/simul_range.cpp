#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include <ctime>
#include <armadillo>
#include <gsl/gsl_rng.h>
#include "geometry.hpp"
#include "utils.hpp"

using namespace std;
using namespace arma;

// Struct to pack the simulation
// parameters all together.
struct Simul_params
{
    // Geometric parameters
    int p;
    int q;
    int dim;

    // Integration parameters
    int L;
    double dt;
    int M;

    // Other parameters
    int iter_therm;
    int iter_simul;
    int gap;
    double g2_i;
    double g2_f;
    double g2_step;

    // HMC mode
    string mode;

    // Validity string
    string valid;
};

// Function to read simulation parameters from stream
bool read_init_stream(istream& in, struct Simul_params& sm);

// Function to check whether parameters have been entered in the
// correct order
bool params_validity(struct Simul_params& sm);


int main(int argc, char** argv)
{
    // ARGUMENT LIST:
    // The first argument is time(NULL) (it will be used as global seed).
    // The second argument is the job rank (it will be used as local seed).
    if(argc != 3)
    {
        cerr << "Error: need to pass global seed and local rank as arguments." << endl;
        return 0;
    }
    // Both arguments are converted to unsigned long.
    unsigned long global_seed = stoul(argv[1]);
    unsigned long local_rank = stoul(argv[2]);
   

    // RNG:
    // Initialize random number generator with global+local
    // (to avoid identical initialization among different jobs).
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(engine, global_seed+local_rank);


    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file yyyymmdd/local_rank/init.tmp
    string foldername = foldername_from_time(global_seed);
    string localpath = foldername + argv[2] + "/";
    string init_filename = localpath + "init.tmp";

    struct Simul_params sm;
    ifstream in_init;
    in_init.open(init_filename);
    
    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + init_filename << endl;
        return 0;
    }

    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " is probably not formatted in the correct way." << endl;
        cerr << "The correct formatting is p:q:dim:L:dt:M:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:mode:" << endl;
        cerr << "Validity string:          " << sm.valid << endl;
        return 0;
    }

    in_init.close();
    
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN MONTE CARLO **********//
    clog << "RNG seed: " << global_seed+local_rank << endl << endl;

    // Open file with precomputed dt values for each g2
    string dt_filename = localpath + "dt.tmp";
    ifstream in_dt;
    in_dt.open(dt_filename);
    
    string prefix = "GEOM";
    double g2;
    while(in_dt >> g2)
    {
        // Create geometry
        Geom24 G(sm.p, sm.q, sm.dim, g2);
        clog << G << endl;


        if(sm.mode == "fix_nosplit")
        {
            // THERMALIZATION
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_nosplit(sm.L, sm.dt, sm.iter_therm, engine, 0.65);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            in_dt >> sm.dt;

            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_fix_nosplit(sm.L, sm.dt, sm.iter_simul, sm.gap, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "fix_split")
        {
            // THERMALIZATION
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_split(sm.L, sm.dt, sm.M, sm.iter_therm, engine, 0.65);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            in_dt >> sm.dt;

            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_fix_split(sm.L, sm.dt, sm.M, sm.iter_simul, sm.gap, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "rand_nosplit")
        {
            // THERMALIZATION
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_nosplit(sm.L, sm.dt, 100, engine, 0.65);
            double dt_min, dt_max;
            in_dt >> dt_min >> dt_max;
            G.HMC_rand_nosplit(sm.L, sm.L, dt_min, dt_max, sm.iter_therm, engine);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_rand_nosplit(sm.L, sm.L, dt_min, dt_max, sm.iter_simul, sm.gap, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "rand_split")
        {
            // THERMALIZATION
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_split(sm.L, sm.dt, sm.M, 100, engine, 0.65);
            double dt_min, dt_max;
            in_dt >> dt_min >> dt_max;
            G.HMC_rand_split(sm.L, sm.L, dt_min, dt_max, sm.M, sm.iter_therm, engine);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_rand_split(sm.L, sm.L, dt_min, dt_max, sm.M, sm.iter_simul, sm.gap, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;

        }
    }

    in_dt.close();
    
    
    //********* END MONTE CARLO **********//


    gsl_rng_free(engine);

    return 0;
}


bool read_init_stream(istream& in, struct Simul_params& sm)
{
    bool success = false;
    if(in)
    {
        string temp;
        
        in >> temp >> sm.p;
        sm.valid += temp;
        in >> temp >> sm.q;
        sm.valid += temp;
        in >> temp >> sm.dim;
        sm.valid += temp;
        in >> temp >> sm.L;
        sm.valid += temp;
        in >> temp >> sm.dt;
        sm.valid += temp;
        in >> temp >> sm.M;
        sm.valid += temp;
        in >> temp >> sm.iter_therm;
        sm.valid += temp;
        in >> temp >> sm.iter_simul;
        sm.valid += temp;
        in >> temp >> sm.gap;
        sm.valid += temp;
        in >> temp >> sm.g2_i;
        sm.valid += temp;
        in >> temp >> sm.g2_f;
        sm.valid += temp;
        in >> temp >> sm.g2_step;
        sm.valid += temp;
        in >> temp >> sm.mode;
        sm.valid += temp;

        success = true;
    }

    return success;
}

bool params_validity(struct Simul_params& sm)
{
    return sm.valid == "p:q:dim:L:dt:M:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:mode:";
}

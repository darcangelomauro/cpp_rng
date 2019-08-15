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
#include "params.hpp"

using namespace std;
using namespace arma;

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
        cerr << "The correct formatting is " << sm.control << endl;
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
            double dt = 0.005;
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_nosplit(sm.L, dt, sm.iter_therm, sm.adj, engine, sm.AR);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            in_dt >> dt;

            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_fix_nosplit(sm.L, dt, sm.iter_simul, sm.gap, sm.adj, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "fix_split")
        {
            // THERMALIZATION
            double dt = 0.005;
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_split(sm.L, dt, sm.M, sm.iter_therm, sm.adj, engine, sm.AR);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            in_dt >> dt;

            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_fix_split(sm.L, dt, sm.M, sm.iter_simul, sm.gap, sm.adj, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "rand_nosplit")
        {
            // THERMALIZATION
            double dt = 0.005;
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_nosplit(sm.L, dt, 100, sm.adj, engine, sm.AR);
            double dt_min, dt_max;
            in_dt >> dt_min >> dt_max;
            G.HMC_rand_nosplit(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.iter_therm, sm.adj, engine);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_rand_nosplit(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.iter_simul, sm.gap, sm.adj, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
        }
        
        else if(sm.mode == "rand_split")
        {
            // THERMALIZATION
            double dt = 0.005;
            clog << "Thermalization start timestamp: " << time(NULL) << endl;
            G.HMC_fix_split(sm.L, dt, sm.M, 100, sm.adj, engine, sm.AR);
            double dt_min, dt_max;
            in_dt >> dt_min >> dt_max;
            G.HMC_rand_split(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.M, sm.iter_therm, sm.adj, engine);
            clog << "Thermalization end timestamp: " << time(NULL) << endl;


            // SIMULATION
            ofstream out_s, out_hl;
            string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            out_s.open(out_filename + "_S.txt");
            out_hl.open(out_filename + "_HL.txt");
            
            clog << "Simulation start timestamp: " << time(NULL) << endl;
            double ar = G.HMC_rand_split(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.M, sm.iter_simul, sm.gap, sm.adj, engine, out_s, out_hl);
            clog << "Simulation end timestamp: " << time(NULL) << endl;

            out_s.close();
            out_hl.close();

            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;

        }
    }

    in_dt.close();
    
    
    //********* END MONTE CARLO **********//


    gsl_rng_free(engine);

    return 0;
}

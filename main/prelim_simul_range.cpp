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
    // The first argument is time(NULL) (it will be used as seed).
    if(argc != 2)
    {
        cerr << "Error: need to pass global seed as arguments." << endl;
        return 0;
    }
    // Convert argument to unsigned long.
    unsigned long global_seed = stoul(argv[1]);
   

    // RNG:
    // Initialize random number generator with global_seed-1
    gsl_rng* engine = gsl_rng_alloc(gsl_rng_ranlxd2);
    gsl_rng_set(engine, global_seed-1);


    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file yyyymmdd/init.txt
    string foldername = foldername_from_time(global_seed);
    string init_filename = foldername + "init.txt";

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



    //********* BEGIN PRELIMINARY MONTE CARLO RUN **********//
    
    // Output file with final dt values for each g2
    string dt_filename = foldername + "dt.txt";
    ofstream out_dt;
    out_dt.open(dt_filename);

    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Create geometry
        Geom24 G(sm.p, sm.q, sm.dim, g2);
        clog << G << endl;


        // PRELIMINARY RUN
        ofstream out_s;
        string out_filename = foldername + "PRELIM/" + filename_from_data(sm.p, sm.q, sm.dim, g2, "PRELIM");
        out_s.open(out_filename + "_S.txt");
        
        clog << "Preliminary run start timestamp: " << time(NULL) << endl;
        
        if(sm.mode == "fix_nosplit")
        {
            // Thermalize
            G.HMC_fix_nosplit(sm.L, sm.dt, sm.iter_therm, engine, 0.65);
            
            // Start run
            double ar = G.HMC_fix_nosplit(sm.L, sm.dt, sm.iter_therm, 1, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << sm.dt << endl;
        }

        else if(sm.mode == "fix_split")
        {
            // Thermalize
            G.HMC_fix_split(sm.L, sm.dt, sm.M, sm.iter_therm, engine, 0.65);
            
            // Start run
            double ar = G.HMC_fix_split(sm.L, sm.dt, sm.M, sm.iter_therm, 1, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << sm.dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << sm.dt << endl;
        }

        else if(sm.mode == "rand_nosplit")
        {
            // Thermalize
            double dt_min = sm.dt;
            G.HMC_fix_nosplit(sm.L, dt_min, sm.iter_therm, engine, 0.9);
            double dt_max = sm.dt;
            G.HMC_fix_nosplit(sm.L, dt_max, sm.iter_therm, engine, 0.4);
            
            // Start run
            double ar = G.HMC_rand_nosplit(sm.L, sm.L, dt_min, dt_max, sm.iter_therm, 1, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << dt_min << " " << dt_max << endl;
        }
        
        else if(sm.mode == "rand_split")
        {
            // Thermalize
            double dt_min = sm.dt;
            G.HMC_fix_split(sm.L, dt_min, sm.M, sm.iter_therm, engine, 0.9);
            double dt_max = sm.dt;
            G.HMC_fix_split(sm.L, dt_max, sm.M, sm.iter_therm, engine, 0.4);
            
            // Start run
            double ar = G.HMC_rand_split(sm.L, sm.L, dt_min, dt_max, sm.M, sm.iter_therm, 1, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << dt_min << " " << dt_max << endl;
        }
        
        
        out_s.close();
        

        // Increment g2
        g2 += sm.g2_step;
    }
    
    //********* END PRELIMINARY MONTE CARLO RUN **********//

    out_dt.close();
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
        in >> sm.mode;
        sm.valid += temp;

        success = true;
    }

    return success;
}

bool params_validity(struct Simul_params& sm)
{
    return sm.valid == "p:q:dim:L:dt:M:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:mode:";
}

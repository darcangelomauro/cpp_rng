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
    double scale;

    // Other parameters
    int iter_therm;
    int iter_simul;
    int gap;
    double g2_i;
    double g2_f;
    double g2_step;

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
        cerr << "The correct formatting is p:q:dim:scale:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:" << endl;
        cerr << "Validity string:          " << sm.valid << endl;
        return 0;
    }

    in_init.close();
    
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN MONTE CARLO **********//
    clog << "RNG seed: " << global_seed+local_rank << endl << endl;

    // Open file with precomputed scale values for each g2
    string scale_filename = localpath + "scale.tmp";
    ifstream in_scale;
    in_scale.open(scale_filename);
    
    string prefix = "GEOM";
    double g2;
    while(in_scale >> g2)
    {
        // Create geometry
        Geom24 G(sm.p, sm.q, sm.dim, g2);
        clog << G << endl;


        // THERMALIZATION
        clog << "Thermalization start timestamp: " << time(NULL) << endl;
        G.MMC(sm.scale, sm.iter_therm, engine, 0.232);
        clog << "Thermalization end timestamp: " << time(NULL) << endl;


        // SIMULATION
        ofstream out_s, out_hl;
        string out_filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
        out_s.open(out_filename + "_S.txt");
        out_hl.open(out_filename + "_HL.txt");
        
        in_scale >> sm.scale;

        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = G.MMC(sm.scale, sm.iter_simul, sm.gap, engine, out_s, out_hl);
        clog << "Simulation end timestamp: " << time(NULL) << endl;

        out_s.close();
        out_hl.close();

        clog << "Metropolis scale: " << sm.scale << endl;
        clog << "Acceptance rate: " << ar << endl;
        clog << endl;

        g2 += sm.g2_step;
    }

    in_scale.close();
    
    
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
        in >> temp >> sm.scale;
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

        success = true;
    }

    return success;
}

bool params_validity(struct Simul_params& sm)
{
    return sm.valid == "p:q:dim:scale:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:";
}

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

    // Other parameters
    int iter_therm;
    int iter_simul;
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


    // FOLDERS:
    // The common folder name for all jobs is of the form yyyymmdd,
    // and should already exist.
    // The output of each job is stored in yyyymmdd/local_rank,
    // which is created here.
    string foldername = foldername_from_time(global_seed);
    string localpath = foldername + argv[2] + "/";
    if(mkdir((localpath).c_str(), 0777))
    {
        cerr << "Error: local folder couldn't be created." << endl;
        cerr << "Possible causes:" << endl;
        cerr << "1) The folder " + foldername + " does not exist." << endl;
        cerr << "2) The folder " + localpath + " already exists." << endl;
        return 0;
    }


    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file yyyymmdd/init.txt
    
    struct Simul_params sm;
    ifstream in_init;
    in_init.open(foldername + "init.txt");
    
    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + foldername + "init.txt" << endl;
        return 0;
    }

    if(!params_validity(sm))
    {
        cerr << "Error: file " + foldername + "init.txt is probably not formatted in the correct way." << endl;
        cerr << "The correct formatting is p:q:dim:L:dt:iter_therm:iter_simul:g2_i:g2_f:g2_step:" << endl;
        cerr << "Validity string:          " << sm.valid << endl;
        return 0;
    }

    in_init.close();
    
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN MONTE CARLO **********//

    double num = (sm.g2_f - sm.g2_i)/sm.g2_step;
    clog << "Starting simulation in g2 range: ";
    clog << "[" << sm.g2_i << " : " << sm.g2_f << "), " << int(num) << " uniformly distributed values" << endl;
    clog << "Local seed: " << global_seed+local_rank << endl;
    clog << endl;
    
    string prefix = "GEOM";
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        Geom24 G(sm.p, sm.q, sm.dim, g2);
       
        clog << G << endl;

        string filename = localpath + filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);


        // THERMALIZATION

        ofstream out_s, out_hl;
        out_s.open(filename + "_S_therm.txt");
        out_hl.open(filename + "_HL_therm.txt");

        clog << "Thermalization start timestamp: " << time(NULL) << endl;
        G.HMC(sm.L, sm.dt, sm.iter_therm, true, engine, out_s, out_hl);
        clog << "Thermalization end timestamp: " << time(NULL) << endl;
        clog << "Integration step: " << sm.dt << endl;
        
        out_s.close();
        out_hl.close();
       
        thermalization_analysis(filename);


        // SIMULATION

        out_s.open(filename + "_S.txt");
        out_hl.open(filename + "_HL.txt");

        clog << "Simulation start timestamp: " << time(NULL) << endl;
        double ar = G.HMC(sm.L, sm.dt, sm.iter_simul, false, engine, out_s, out_hl);
        clog << "Simulation end timestamp: " << time(NULL) << endl;

        out_s.close();
        out_hl.close();


        clog << "Acceptance rate: " << ar << endl;
        clog << endl;

        g2 += sm.g2_step;
    }
    
    
    //********* END MONTE CARLO **********//


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
        in >> temp >> sm.iter_therm;
        sm.valid += temp;
        in >> temp >> sm.iter_simul;
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
    return sm.valid == "p:q:dim:L:dt:iter_therm:iter_simul:g2_i:g2_f:g2_step:";
}

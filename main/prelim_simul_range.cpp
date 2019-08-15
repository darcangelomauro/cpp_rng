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
        cerr << "The correct formatting is " << sm.control << endl;
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
            double dt = 0.005;
            G.HMC_fix_nosplit(sm.L, dt, sm.iter_therm, sm.adj, engine, sm.AR);
            
            // Start run
            double ar = G.HMC_fix_nosplit(sm.L, dt, sm.iter_therm, 1, sm.adj, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << dt << endl;
        }

        else if(sm.mode == "fix_split")
        {
            // Thermalize
            double dt = 0.005;
            G.HMC_fix_split(sm.L, dt, sm.M, sm.iter_therm, sm.adj, engine, sm.AR);
            
            // Start run
            double ar = G.HMC_fix_split(sm.L, dt, sm.M, sm.iter_therm, 1, sm.adj, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_dt << g2 << " " << dt << endl;
        }

        else if(sm.mode == "rand_nosplit")
        {
            // Thermalize
            double dt_min = 0.005;
            G.HMC_fix_nosplit(sm.L, dt_min, sm.iter_therm, sm.adj, engine, sm.AR+sm.dAR);
            double dt_max = 0.005;
            G.HMC_fix_nosplit(sm.L, dt_max, sm.iter_therm, sm.adj, engine, sm.AR-sm.dAR);
            
            // Start run
            double ar = G.HMC_rand_nosplit(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.iter_therm, 1, sm.adj, engine, out_s);
            
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
            double dt_min = 0.005;
            G.HMC_fix_split(sm.L, dt_min, sm.M, sm.iter_therm, sm.adj, engine, sm.AR+sm.dAR);
            double dt_max = 0.005;
            G.HMC_fix_split(sm.L, dt_max, sm.M, sm.iter_therm, sm.adj, engine, sm.AR-sm.dAR);
            
            // Start run
            double ar = G.HMC_rand_split(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.M, sm.iter_therm, 1, sm.adj, engine, out_s);
            
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

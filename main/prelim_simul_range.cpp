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
        return 1;
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
        return 1;
    }

    clog << "File " + init_filename + " contains the following parameters:" << endl;
    clog << sm.control << endl;

    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
        return 1;
    }

    in_init.close();
    
    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN PRELIMINARY MONTE CARLO RUN **********//
    
    // Output file with final dt or scale values for each g2
    string duav_filename = foldername + "duav.txt";
    ofstream out_duav;
    out_duav.open(duav_filename);

    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Create geometry
        Geom24 G(sm.p, sm.q, sm.dim, g2);
        G.shuffle(engine);
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
            out_duav << g2 << " " << dt << endl;
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
            out_duav << g2 << " " << dt << endl;
        }

        else if(sm.mode == "rand_nosplit")
        {
            // Thermalize
            double dt_min = 0.005;
            G.HMC_fix_nosplit(sm.L, dt_min, sm.iter_therm, sm.adj, engine, sm.AR+sm.dAR);
            double dt_max = 0.005;
            G.HMC_fix_nosplit(sm.L, dt_max, sm.iter_therm, sm.adj, engine, sm.AR-sm.dAR);

            if(dt_min > dt_max)
            {
                double temp = dt_min;
                dt_min = dt_max;
                dt_max = temp;
            }
            
            // Start run
            double ar = G.HMC_rand_nosplit(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.iter_therm, 1, sm.adj, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_duav << g2 << " " << dt_min << " " << dt_max << endl;
        }
        
        else if(sm.mode == "rand_split")
        {
            // Thermalize
            double dt_min = 0.005;
            G.HMC_fix_split(sm.L, dt_min, sm.M, sm.iter_therm, sm.adj, engine, sm.AR+sm.dAR);
            double dt_max = 0.005;
            G.HMC_fix_split(sm.L, dt_max, sm.M, sm.iter_therm, sm.adj, engine, sm.AR-sm.dAR);
            
            if(dt_min > dt_max)
            {
                double temp = dt_min;
                dt_min = dt_max;
                dt_max = temp;
            }
            
            // Start run
            double ar = G.HMC_rand_split(sm.L-sm.dL, sm.L+sm.dL, dt_min, dt_max, sm.M, sm.iter_therm, 1, sm.adj, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Integration step: " << dt_min << " " << dt_max << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of dt
            out_duav << g2 << " " << dt_min << " " << dt_max << endl;
        }
        
        else if(sm.mode == "mmc")
        {
            // Thermalize
            double scale = 0.005;
            G.MMC(scale, sm.iter_therm, sm.adj, engine, sm.AR);
            
            // Start run
            double ar = G.MMC(scale, sm.iter_therm, 1, sm.adj, engine, out_s);
            
            // Output log
            clog << "Preliminary run end timestamp: " << time(NULL) << endl;
            clog << "Metropolis scale: " << scale << endl;
            clog << "Acceptance rate: " << ar << endl;
            clog << endl;
            
            // Write final value of scale
            out_duav << g2 << " " << scale << endl;
        }
        
        out_s.close();
        

        // Increment g2
        g2 += sm.g2_step;
    }
    
    //********* END PRELIMINARY MONTE CARLO RUN **********//

    out_duav.close();
    gsl_rng_free(engine);

    return 0;
}

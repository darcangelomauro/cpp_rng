#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "geometry.hpp"
#include "utils.hpp"
#include "params.hpp"
#include "statistics.hpp"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    // Check arguments
    if(argc < 5)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of the folder containing the data" << endl;
        cerr << "2) First index of the jobs array" << endl;
        cerr << "3) Number of jobs in the array" << endl;
        cerr << "4) Name of the observable" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    int fst_jarr = stoi(argv[2]);
    int num_jarr = stoi(argv[3]);
    string name = argv[4]; 



    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string init_filename = path + "/init.txt";

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



    //********* BEGIN ANALYSIS **********//
    

    // Open output file 
    string out_filename = path + "/observables/" + name + ".txt";
    ofstream out_obs;
    out_obs.open(out_filename);

    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Create vector of uncorrelated samples
        vec samples(num_jarr);

        // Cycle on jobs in the array
        for(int i=0; i<num_jarr; ++i)
        {
            // Open data files
            string array_path = path + "/" + to_string(i+fst_jarr);
            string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_s;
            in_s.open(array_path + "/" + filename + "_S.txt");

            // Create vector of correlated samples
            int length = n_meas(sm.iter_simul, sm.gap);
            vec vec_corr(length);

            // Cycle on samples
            for(int j=0; j<length; ++j) 
            {
                double S2, S4;
                in_s >> S2 >> S4;

                // ***** COMPUTE OBSERVABLE HERE *****
                double temp = 2*g2*S2 + 4*S4;
                // ***** THAT'S IT, YOU'RE DONE *****

                vec_corr(j) = temp;
            }
            in_s.close();

            // Initialize i-th element of vector of uncorrelated samples with mean of job #i
            samples(i) = mean(vec_corr);
        }

        // Find theoretical value
        Geom24 G(sm.p, sm.q, 1, 1);
        double c = G.get_nHL()*sm.dim*sm.dim;

        if(num_jarr > 1)
        {
            // Output mean and error of observable
            double avg = 0;
            double var = 0;
            jackknife(samples, avg, var, my_mean);
            double err = sqrt(var);
            out_obs << g2 << " " << avg << " " << err << " " << c << endl;
        }
        else
            out_obs << g2 << " " << samples(0) << " " << "---" << " " << c << endl;


        g2 += sm.g2_step;
    }

    out_obs.close();

    //********* END ANALYSIS **********//

    return 0;
}

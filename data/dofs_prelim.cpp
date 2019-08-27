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
    if(argc < 2)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of the folder containing the data" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "PRELIM";
    string path = argv[1];



    //********* BEGIN PARAMETER INITIALIZATION **********//
    
    // Read simulation parameters from file path/init.txt
    string init_filename = path + "/init.txt";

    struct Simul_params sm;
    ifstream in_init;
    in_init.open(init_filename);

    if(!read_init_stream(in_init, sm))
    {
        cerr << "Error: couldn't read file " + init_filename << endl;
        return 1;
    }

    cout << "File " + init_filename + " contains the following parameters:" << endl;
    cout << sm.control << endl;

    if(!params_validity(sm))
    {
        cerr << "Error: file " + init_filename + " does not contain the necessary parameters." << endl;
        return 1;
    }

    in_init.close();

    //********* END PARAMETER INITIALIZATION **********//



    //********* BEGIN ANALYSIS **********//
    

    // Open output file 
    string out_filename = path + "/observables/dofs_prelim.txt";
    ofstream out_obs;
    out_obs.open(out_filename);

    if(!out_obs)
    {
        cerr << "Error: file " + out_filename + " could not be opened." << endl;
        return 1;
    }

    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Open data files
        string array_path = path + "/PRELIM";
        string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
        ifstream in_s;
        in_s.open(array_path + "/" + filename + "_S.txt");

        // Create vector of correlated samples
        int length = sm.iter_therm;
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

        // Find theoretical value
        Geom24 G(sm.p, sm.q, 1, 1);
        double c = G.get_nHL()*sm.dim*sm.dim - G.get_nL();

        // Output mean and error of observable
        double avg = 0;
        double var = 0;
        jackknife(vec_corr, avg, var, my_mean);
        double err = sqrt(var);
        out_obs << g2 << " " << avg << " " << err << " " << c << endl;

        g2 += sm.g2_step;
    }

    out_obs.close();

    //********* END ANALYSIS **********//

    return 0;
}

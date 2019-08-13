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
    


    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Cycle on jobs in the array
        for(int i=0; i<num_jarr; ++i)
        {
            // Open data files
            string array_path = path + "/" + to_string(i+fst_jarr);
            string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_s, in_hl;
            in_s.open(array_path + "/" + filename + "_S.txt");
            in_hl.open(array_path + "/" + filename + "_HL.txt");

            // Create vector of correlated samples
            int length = n_meas(sm.iter_simul, sm.gap);

            // Open output file 
            string coord_filename = path + "/observables/" + to_string(i) + "_" + to_string(g2) + "_coord.txt";
            ofstream out_coord;
            out_coord.open(coord_filename);
            
            // Cycle on samples
            for(int j=0; j<length; ++j) 
            {
                Geom24 G(sm.p, sm.q, sm.dim, g2);
                double S2, S4;
                in_s >> S2 >> S4;
                G.read_mat(in_hl);

                // ***** COMPUTE OBSERVABLE HERE *****
                double temp1 = trace(G.get_mat(0)).real();
                double temp2 = trace(G.get_mat(1)).real();
                // ***** THAT'S IT, YOU'RE DONE *****

                out_coord << temp1 << " " << temp2 << endl;
            }
            in_s.close();
            in_hl.close();
            out_coord.close();

        }


        g2 += sm.g2_step;
    }

    //********* END ANALYSIS **********//

    return 0;
}

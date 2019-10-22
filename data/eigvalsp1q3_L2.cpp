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
    if(argc < 4)
    {
        cerr << "Need to pass:" << endl;
        cerr << "1) Name of the folder containing the data" << endl;
        cerr << "2) First index of the jobs array" << endl;
        cerr << "3) Number of jobs in the array" << endl;
        return 1;
    }

    // Some declarations for later
    string prefix = "GEOM";
    string path = argv[1];
    int fst_jarr = stoi(argv[2]);
    int num_jarr = stoi(argv[3]);



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
    


    // Cycle on g2 values
    double g2 = sm.g2_i;
    while(g2 < sm.g2_f)
    {
        // Print value of g2 being precessed
        clog << "g2: " << g2 << endl;

        // Cycle on jobs in the array
        for(int i=0; i<num_jarr; ++i)
        {
            string array_path = path + "/" + to_string(i+fst_jarr);
            Geom24 temp(sm.p, sm.q, 1, 1);
            
            // Compute number of measurements
            int length = n_meas(sm.iter_simul, sm.gap);
            
            // eigenvalues matrix
            mat M(sm.dim, length);
            
            // Open output files
            string out_filename = array_path + "/" + filename_from_data(sm.p, sm.q, sm.dim, g2, "L2") + ".txt";
            ofstream out_obs;
            out_obs.open(out_filename);
            
            if(!out_obs)
            {
                cerr << "Error: file " + out_filename + " could not be opened." << endl;
                return 1;
            }

            // Open data files
            string filename = filename_from_data(sm.p, sm.q, sm.dim, g2, prefix);
            ifstream in_s, in_hl;
            in_s.open(array_path + "/" + filename + "_S.txt");
            in_hl.open(array_path + "/" + filename + "_HL.txt");


            
            // Cycle on measurements
            for(int j=0; j<length; ++j) 
            {
                Geom24 G(sm.p, sm.q, sm.dim, g2);
                double S2, S4;
                in_s >> S2 >> S4;
                G.read_mat(in_hl);


                // ***** COMPUTE OBSERVABLE HERE *****
                vec eigval = eig_sym(G.get_mat(3));
                for(int l=0; l<G.get_dim(); ++l)
                {
                    M(l,j) = eigval(l);
                }
            }
            in_s.close();
            in_hl.close();



            
            for(uword a=0; a<M.n_rows; ++a)
            {
                for(uword b=0; b<M.n_cols; ++b)
                    out_obs << M(a,b) << " ";
                
                out_obs << endl;
            }
            out_obs.close();

        }

        g2 += sm.g2_step;
    }

    //********* END ANALYSIS **********//

    return 0;
}

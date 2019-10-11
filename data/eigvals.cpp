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
            
            // Open output files
            string* out_filename = new string [temp.get_nHL() + 1];
            ofstream* out_obs = new ofstream [temp.get_nHL() + 1];
            mat* M = new mat [temp.get_nHL() + 1];
            for(int j=0; j<temp.get_nHL()+1; ++j)
            {
                string mat_name;
                if(j < temp.get_nH())
                {
                    mat_name = "EVALSH" + to_string(j);
                    M[j] = mat(sm.dim, length); 
                }
                else if(j == temp.get_nHL())
                {
                    mat_name = "EVALSD";
                    M[j] = mat(sm.dim*sm.dim*temp.get_dim_omega(), length); 
                }
                else
                {
                    mat_name = "EVALSL" + to_string(j-temp.get_nH());
                    M[j] = mat(sm.dim, length); 
                }

                out_filename[j] = array_path + "/" + filename_from_data(sm.p, sm.q, sm.dim, g2, mat_name) + ".txt";
                out_obs[j].open(out_filename[j]);
            
                if(!out_obs)
                {
                    cerr << "Error: file " + out_filename[j] + " could not be opened." << endl;
                    return 1;
                }
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
                for(int k=0; k<G.get_nHL(); ++k)
                {
                    vec eigval = eig_sym(G.get_mat(k));
                    for(int l=0; l<G.get_dim(); ++l)
                    {
                        M[k](l,j) = eigval(l);
                    }
                }

                cx_mat D = G.build_dirac();
                vec eigval = eig_sym(D);
                for(uword l=0; l<eigval.n_elem; ++l)
                {
                    M[G.get_nHL()](l,j) = eigval(l);
                }
            }
            in_s.close();
            in_hl.close();



            
            for(int j=0; j<temp.get_nHL()+1; ++j)
            {
                for(uword a=0; a<M[j].n_rows; ++a)
                {
                    for(uword b=0; b<M[j].n_cols; ++b)
                        out_obs[j] << M[j](a,b) << " ";
                    
                    out_obs[j] << endl;
                }
                out_obs[j].close();
            }

            delete [] out_filename; 
            delete [] out_obs; 
            delete [] M; 
        }

        g2 += sm.g2_step;
    }

    //********* END ANALYSIS **********//

    return 0;
}

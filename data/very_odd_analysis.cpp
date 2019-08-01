#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
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
    int gap;
    double g2_i;
    double g2_f;
    double g2_step;

    // Validity string
    string valid;
};

// Function to read simulation parameters from stream
bool read_init_stream(istream&, struct Simul_params&);

// Function to check whether parameters have been entered in the
// correct order
bool params_validity(struct Simul_params&);

// Jack knife mean and variance estimate of an arbitrary function f
void jackknife(const vec&, double&, double&, double f(const vec&));

// Some stat functions
double my_mean(const vec&);
double my_var(const vec&);
double my_sus(const vec&);

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
        cerr << "The correct formatting is p:q:dim:L:dt:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:" << endl;
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
    return sm.valid == "p:q:dim:L:dt:iter_therm:iter_simul:gap:g2_i:g2_f:g2_step:";
}


void jackknife(const vec& vec_uncorr, double& avg, double& var, double f(const vec&))
{
    // Find vector size
    int size = vec_uncorr.n_elem;

    // Create vector of delete-1 clusters
    vec vec_del1(size);
    
    // Calculate delete-1 clusters
    for(int i=0; i<size; ++i)
    {
        // Create i-th cluster
        vec ith_cluster(size-1);
        for(int j=0; j<size; ++j)
        {
            // Copy element in the same position if j < i
            if(j < i)
                ith_cluster(j) = vec_uncorr(j);
            // Copy element one position behind if j > i
            else if(j > i)
                ith_cluster(j-1) = vec_uncorr(j);

        }

        // Put the return value of f into vec_del1 in i-th position
        vec_del1(i) = f(ith_cluster);
    }

    // Calclulate unbiased mean
    avg = mean(vec_del1);

    // Calculate unbiased variance
    var = 0;
    for(int i=0; i<size; ++i)
        var += pow(vec_del1(i) - avg, 2);
    var *= (double)(size-1)/size;
}

double my_mean(const vec& vec_uncorr)
{
    return mean(vec_uncorr);
}

double my_var(const vec& vec_uncorr)
{
    return var(vec_uncorr, 1);
}

double my_sus(const vec& vec_uncorr)
{
    return var(vec_uncorr, 1)/mean(vec_uncorr);
}

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
    

    // Open output file 
    string out_filename = path + "/observables/dofs_prelim.txt";
    ofstream out_obs;
    out_obs.open(out_filename);

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
        double c = G.get_nHL()*sm.dim*sm.dim;

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

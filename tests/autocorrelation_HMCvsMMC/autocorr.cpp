#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "geometry.hpp"

using namespace std;
using namespace arma;


int main(int argc, char** argv)
{
    if(argc < 2)
    {
        cout << "Too few arguments. Needs dataset filename. Quitting" << endl;
        return 1;
    }
    string argv1(argv[1]);

    // open autocorrelation file
    ofstream out_auto;
    out_auto.open("data/" + argv1 + "_autocorr.txt");
    
    // read observable from file
    ifstream in_obs;
    in_obs.open("data/" + argv1 + ".txt");
    if(in_obs.is_open())
    {
        // fval will store the content of the file
        vector<double> fval;

        // read from file
        double temp;
        while(in_obs >> temp)
            fval.push_back(temp);

        // find vector size (i.e. number of measurements)
        vector<double>::size_type Nsw;
        Nsw = fval.size();

        // compute autocorrelation and store in autval
        vector<double> autval;
        for(unsigned i=0; i<Nsw; ++i)
        {
            double term1=0;
            double term2=0;
            double term3=0;
            
            double length = Nsw-i;
            for(unsigned j=0; j<length; ++j)
            {
                term1 += fval[j]*fval[i+j];
                term2 += fval[j];
                term3 += fval[i+j];
            }


            autval.push_back(term1/length - term2*term3/(length*length));
        }

        // print on file
        vector<double>::const_iterator begin = autval.begin();
        vector<double>::const_iterator end = autval.end();
        for(vector<double>::const_iterator iter = begin; iter != end; ++iter)
            out_auto << (*iter)/(*begin) << endl;

    }
    else
    {
        cout << "File data/" << argv1 << ".txt not found. Quitting" << endl;
        return 1;
    }
    in_obs.close();
    out_auto.close();

    return 0;
}

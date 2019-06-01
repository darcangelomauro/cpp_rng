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
    if(argc < 3)
    {
        cout << "Too few arguments. Needs dataset filename and either h (for hmc) or m (for mmc). Quitting" << endl;
        return 1;
    }

    // open data file
    string argv1(argv[1]);
    string horm(argv[2]);
    string filename_data = "data/" + argv1 + "_data.txt";
    ifstream in_data;
    in_data.open(filename_data);
    
    if(in_data.is_open())
    {
        // open observable file
        ofstream out_f;
        out_f.open("data/" + argv1 + "_F_" + horm + ".txt");
        
        // read data from data file
        int p, q, dim, Nsw;
        double g;
        in_data >> p >> q >> dim >> g >> Nsw;

        // read the matrices from file
        string filename_hl = "data/" + argv1 + "_HL_" + horm + ".txt";
        ifstream in_hl;
        in_hl.open(filename_hl);
        if(in_hl.is_open())
        {
            for(int i=0; i<Nsw; ++i)
            {
                Geom24 G(p, q, dim, g);
                G.read_mat(in_hl);

                // compute observable
                double f = 0;
                double den = 0;
                for(int j=0; j<p; ++j)
                {
                    double trH = trace(G.get_mat(j)).real();
                    double trH2 = trace(G.get_mat(j)*G.get_mat(j)).real();
                    f += trH*trH;
                    den += trH2;
                }
                
                // output observable
                out_f << f/(dim*den) << endl;
            }
        }
        else
        {
            cout << "File " << filename_hl << " not found. Quitting" << endl;
            return 1;
        }
        in_hl.close();
        out_f.close();

    }
    else
    {
        cout << "File " << filename_data << " not found. Quitting" << endl;
        return 1;
    }
    in_data.close();



    return 0;
}

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
    string filename_data = "data/" + argv1 + "_data.txt";
    ifstream in_data;
    in_data.open(filename_data);
    
    if(in_data.is_open())
    {
        ofstream out_s;
        out_s.open("data/" + argv1 + "_S_varG.txt");
        ofstream out_f;
        out_f.open("data/" + argv1 + "_F_varG.txt");
        
        int idx;
        int p, q, dim;
        double g;
        while(in_data >> idx >> p >> q >> dim >> g)
        {
            // action related operations
            string filename_s = "data/" + argv1 + "_simS_" + to_string(idx) + ".txt";
            ifstream in_s;
            in_s.open(filename_s);
            
            vector<double>::size_type Nsw;
            if(in_s.is_open())
            {
                vector<double> sval;
                double s, k, h;
                while(in_s >> s >> k >> h)
                    sval.push_back(s);
                in_s.close();

                Nsw = sval.size();
                vector<double>::const_iterator begin(sval.begin());
                vector<double>::const_iterator end(sval.end());

                out_s << g << " " << accumulate(begin+Nsw/2, end, 0.0)/(Nsw/2) << endl;
            }
            else
            {
                cout << "File " << filename_s << " not found. Quitting" << endl;
                return 1;
            }
            in_s.close();

            // matrix related operations
            string filename_hl = "data/" + argv1 + "_simHL_" + to_string(idx) + ".txt";
            ifstream in_hl;
            in_hl.open(filename_hl);
            if(in_hl.is_open())
            {
                vector<double> fval;
                for(unsigned i=0; i<Nsw; ++i)
                {
                    Geom24 G(p, q, dim, g);
                    G.read_mat(in_hl);

                    double f = 0;
                    double den = 0;
                    for(int j=0; j<p; ++j)
                    {
                        double trH = trace(G.get_mat(j)).real();
                        double trH2 = trace(G.get_mat(j)*G.get_mat(j)).real();
                        f += trH*trH;
                        den += trH2;
                    }
                    fval.push_back(f/(dim*den));
                
                }
                vector<double>::const_iterator begin(fval.begin());
                vector<double>::const_iterator end(fval.end());

                out_f << g << " " << accumulate(begin+Nsw/2, end, 0.0)/(Nsw/2) << endl;
            }
            else
            {
                cout << "File " << filename_hl << " not found. Quitting" << endl;
                return 1;
            }
            in_hl.close();
        }
        out_s.close();
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

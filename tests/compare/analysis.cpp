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
    string filename_data = argv1 + "_data.txt";
    ifstream in_data;
    in_data.open(filename_data);
    
    if(in_data.is_open())
    {
        ofstream out;
        out.open(argv1 + "_varG.txt");
        
        int idx;
        int p, q, dim;
        double g;
        while(in_data >> idx >> p >> q >> dim >> g)
        {
            // action related operations
            string filename_s = argv1 + "_simS_" + to_string(idx) + ".txt";
            ifstream in_s;
            in_s.open(filename_s);
            if(in_s.is_open())
            {
                vector<double> sval;
                double s, k, h;
                while(in_s >> s >> k >> h)
                    sval.push_back(s);
                in_s.close();

                vector<double>::size_type size = sval.size();
                vector<double>::const_iterator begin(sval.begin());
                vector<double>::const_iterator end(sval.end());

                out << g << " " << accumulate(begin+size/2, end, 0.0)/(size/2) << endl;
            }
            else
            {
                cout << "File " << filename_s << " not found. Quitting" << endl;
                return 1;
            }
        }

        out.close();
    }
    else
    {
        cout << "File " << filename_data << " not found. Quitting" << endl;
        return 1;
    }



    return 0;
}

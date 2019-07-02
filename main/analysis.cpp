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


int main(int argc, char** argv)
{
    if(argc < 7)
    {
        cout << "Too few arguments. Needs p, q, dim, start g2, end g2, step g2" << endl;
        return 1;
    }

    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int dim = stoi(argv[3]);
    string basename = basename_from_data(p, q, dim);
    
    double g2 = stod(argv[4]);
    double g2_end = stod(argv[5]);
    double g2_step = stod(argv[6]);
    
    ofstream out_s;
    out_s.open("data/" + basename + "_S_range.txt");
    ofstream out_f;
    out_f.open("data/" + basename + "_F_range.txt");
    ofstream out_comm;
    out_comm.open("data/" + basename + "_COMM_range.txt");

    
    while(g2 <= g2_end)
    {
        ifstream in_s, in_hl;
        string filename = filename_from_data(p, q, dim, g2);
        in_s.open("data/" + filename + "_S.txt");
        in_hl.open("data/" + filename + "_HL.txt");

        if(in_s.is_open() && in_hl.is_open())
        {
            cout << "Processing " + filename << endl;

            // action 
            vector<double>::size_type Nsw;
            vector<double> sval;
            double s2, s4;
            while(in_s >> s2 >> s4)
                sval.push_back(g2*s2+s4);
            in_s.close();

            Nsw = sval.size();
            vector<double>::const_iterator begin(sval.begin());
            vector<double>::const_iterator end(sval.end());

            out_s << g2 << " " << accumulate(begin, end, 0.0)/Nsw << endl;

            // matrix related operations
            vector<double> fval;
            vector<double> commval;
            for(unsigned i=0; i<Nsw; ++i)
            {
                Geom24 G(p, q, dim, g2);
                G.read_mat(in_hl);
                
                // F
                double f = 0;
                double den = 0;
                for(int j=0; j<G.get_nH(); ++j)
                {
                    double trH = trace(G.get_mat(j)).real();
                    double trH2 = trace(G.get_mat(j)*G.get_mat(j)).real();
                    f += trH*trH;
                    den += trH2;
                }
                fval.push_back(f/(dim*den));

                // COMM
                cx_mat L1L1 = G.get_mat(G.get_nH())*G.get_mat(G.get_nH());
                cx_mat L2L2 = G.get_mat(G.get_nH()+1)*G.get_mat(G.get_nH()+1);
                cx_mat L3L3 = G.get_mat(G.get_nH()+2)*G.get_mat(G.get_nH()+2);
                cx_mat L1L2 = G.get_mat(G.get_nH())*G.get_mat(G.get_nH()+1);
                cx_mat L2L1 = G.get_mat(G.get_nH()+1)*G.get_mat(G.get_nH());
                
                double num = trace( (L1L2-L2L1)*G.get_mat(G.get_nH()+2) ).imag();
                num *= num;

                den = 2*trace(L1L1*L2L2 - L1L2*L1L2).real()*trace(L3L3).real();
                den *= den;

                commval.push_back(-num/den);
            }
            in_hl.close();
            
            begin = fval.begin();
            end = fval.end();
            out_f << g2 << " " << accumulate(begin, end, 0.0)/Nsw << endl;
            
            begin = commval.begin();
            end = commval.end();
            out_comm << g2 << " " << accumulate(begin, end, 0.0)/Nsw << endl;
        }

        g2 += g2_step;
    }
    
    out_s.close();
    out_f.close();
    out_comm.close();


    return 0;
}

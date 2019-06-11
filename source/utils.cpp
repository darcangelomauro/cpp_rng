#include <sstream>
#include <iostream>
#include <algorithm>
#include "utils.hpp"

using namespace std;


string filename_from_data(const int& p, const int& q, const int& dim, const double& g2)
{
    string sp = to_string(p);
    string sq = to_string(q);
    string sdim = to_string(dim);
    
    ostringstream osg2;
    osg2 << g2;
    string sg2 = osg2.str();
    replace(sg2.begin(), sg2.end(), '.', 'd');
    
    return "p" + sp + "q" + sq + "dim" + sdim + "g" + sg2;
}

void data_from_filename(const string& s, int& p, int& q, int& dim, double& g2)
{
    string s_p = s.substr(s.find("p")+1, s.find("q")-1);
    string s_q = s.substr(s.find("q")+1, s.find("d")-s.find("q")-1);
    string s_dim = s.substr(s.find("m")+1, s.find("g")-s.find("m")-1);
    string s_g = s.substr(s.find("g")+1, s.find(".")-s.find("g")-1);
    replace(s_g.begin(), s_g.end(), 'd', '.');

    p = stoi(s_p);
    q = stoi(s_q);
    dim = stoi(s_dim);
    g2 = stod(s_g);
}



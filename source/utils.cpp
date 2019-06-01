#include <sstream>
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



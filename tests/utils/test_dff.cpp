#include <string>
#include <iostream>
#include "geometry.hpp"
#include "utils.hpp"

using namespace std;


int main()
{
    int p, q, dim;
    double g2;

    string s = "awfasfabbmvkflrdoeGEOMp1q3dim64g-3d767ddaaaddda.txt";
    string prefix = "GEOM";
    
    data_from_filename(s, p, q, dim, g2, prefix);

    cout << p << " " << q << " " << dim << " " << g2 << endl;

    return 0;
}

#include <string>
#include <iostream>
#include "utils.hpp"

using namespace std;


int main()
{
    int p, q, dim;
    double g2;

    string s = "p1q3dim64g-3d767.txt";
    
    data_from_filename(s, p, q, dim, g2);

    cout << p << " " << q << " " << dim << " " << g2 << endl;

    return 0;
}

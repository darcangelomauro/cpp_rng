// g++ main_cliff.cpp geometry.cpp clifford.cpp -o main_cliff -O2 -std=c++11 -larmadillo
#include <iostream>
#include <cstdlib>
#include "geometry.hpp"
#include "clifford.hpp"
#include <ctime>

using namespace std;


int main()
{
    int p, q;

    cin >> p >> q;

    Cliff C(p, q);

    cout << C << endl;

    cout << "gammas:" << endl;
    for(int i=0; i<C.get_p()+C.get_q(); ++i)
        cout << C.get_gamma(i) << endl;
    cout << "chirality:" << endl;
    cout << C.get_chiral() << endl;
    

    return 0;
}

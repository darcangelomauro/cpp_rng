// g++ main_omega.cpp geometry.cpp clifford.cpp -o main_omega -O2 -std=c++11 -larmadillo
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <armadillo>
#include "geometry.hpp"

using namespace std;
using namespace arma;


int main()
{
    Geom24 G(cin);

    int dim = G.get_dim_omega();

    cout << "nH: " <<G.get_nH() << endl << "nL: " << G.get_nL() << endl;
    
    cout << "omegas:" << endl;
    cout.precision(1);
    for(int i=0; i<G.get_nHL(); ++i)
    {
        cx_mat M = G.get_omega(i);
        
        cout << "omega[" << i << "]:" << endl;
        if(M.is_hermitian())
        {
            cout << "hermitian" << endl;
            if(approx_equal(M*M, cx_mat(dim, dim, fill::eye), "absdiff", 0.0000001))
                cout << "squares to 1" << endl;
        }
        else
            cout << "not hermitian" << endl;
        M.raw_print();
        cout << endl;
    }


    return 0;
}

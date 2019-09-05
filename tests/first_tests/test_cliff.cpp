#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "geometry.hpp"
#include "clifford.hpp"
#include <ctime>

using namespace std;
using namespace arma;


int main()
{
    int p, q;

    cin >> p >> q;

    Cliff C(p, q);
    int dim = C.get_dim_gamma();

    cout << "gammas:" << endl;
    cout.precision(1);
    for(int i=0; i<C.get_p()+C.get_q(); ++i)
    {
        cx_mat M = C.get_gamma(i);
        
        cout << "gamma[" << i << "]:" << endl;
        if(M.is_hermitian())
        {
            cout << "hermitian" << endl;
            if(approx_equal(M*M, cx_mat(dim, dim, fill::eye), "absdiff", 0.0000001))
                cout << "squares to 1" << endl;
        }
        else if((sp_cx_mat(M.t()+M)).n_nonzero == 0)
        {
            cout << "antihermitian" << endl;
            if(approx_equal(M*M, -cx_mat(dim, dim, fill::eye), "absdiff", 0.0000001))
                cout << "squares to -1" << endl;
        }
        else
            cout << "not hermitian nor antihermitian" << endl;
        M.raw_print();
        cout << endl;
    }

    cout << "chirality (derived):" << endl;
    C.get_chiral().raw_print();
    cout << endl;

    cout << "chirality (computed):" << endl;
    cx_mat chiral(dim, dim, fill::eye);
    for(int i=0; i<C.get_p()+C.get_q(); ++i)
        chiral *= C.get_gamma(i);
    int s = (C.get_q() - C.get_p()) % 8;
    int exp = s*(s+1)/2;
    chiral *= pow(cx_double(0,1), exp);
    chiral.raw_print();
    cout << endl;

    cout << "non-vanishing anticommutators:" << endl;
    for(int i=0; i<C.get_p()+C.get_q(); ++i)
    {
        for(int j=0; j<C.get_p()+C.get_q(); ++j)
        {
            cx_mat M1 = C.get_gamma(i);
            cx_mat M2 = C.get_gamma(j);

            cx_mat A = M1*M2 + M2*M1;
            if(!approx_equal(A, cx_mat(dim , dim, fill::zeros), "absdiff", 0.00000001))
                    cout << "{" << i << "," << j << "}" << endl;
        }
    }


    
    

    return 0;
}

// g++ main.cpp geometry.cpp -o main -O2 -std=c++11 -larmadillo
#include <iostream>
#include <cstdlib>
#include "geometry.hpp"
#include <ctime>

using namespace std;


int main()
{
    arma::arma_rng::set_seed(time(NULL));

    Geom24 c1(cin);

    cout << "c1nH: " <<c1.get_nH() << " c1nL: " << c1.get_nL() << endl;

    Geom24 c2(2, 1, 3, 3.453);

    cout << c1 << endl;
    cout << c2 << endl;

    Geom24 c3 = c1;

    c1 = c2;

    cout << c1 << endl;
    cout << c2 << endl;
    cout << c3 << endl;

    cout << c1.MAT(0) << endl;
    cout << "S = " << c1.calculate_S() << endl;
    c1.shuffle();
    cout << c1.MAT(0) << endl;
    cout << "S = " << c1.calculate_S() << endl;


    /*
    for(int i=0; i<10000000; i++)
    {
        Geom24 c4(2, 1, 100, 6.657);
    }
    */
    

    return 0;
}

// g++ main.cpp -o main -O2 -std=c++11
#include <iostream>
#include <cstdlib>
#include "geometry.hpp"

using namespace std;


int main()
{
    Core c1(cin);

    Core c2(2, 1, 56);

    cout << c1 << endl;
    cout << c2 << endl;

    Core c3 = c1;

    c1 = c2;

    cout << c1 << endl;
    cout << c2 << endl;
    cout << c3 << endl;

    

    return 0;
}

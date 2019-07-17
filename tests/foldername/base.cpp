#include <ctime>
#include <iostream>
#include <string>
#include "utils.hpp"

using namespace std;

int main()
{
	time_t t = time(NULL);
    string t_s = to_string(t);
    unsigned long t_ul = stoul(t_s);

    cout << foldername_from_time(t) << endl;
    cout << foldername_from_time(t_ul) << endl;
}

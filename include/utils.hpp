#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

std::string filename_from_data(const int&, const int&, const int&, const double&);
std::string basename_from_data(const int&, const int&, const int&);
void data_from_filename(const std::string&, int&, int&, int&, double&);

#endif


#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

std::string filename_from_data(const int&, const int&, const int&, const double&, const std::string&);
std::string basename_from_data(const int&, const int&, const int&, const std::string&);
void data_from_filename(const std::string&, int&, int&, int&, double&, const std::string&);
std::string foldername_from_time(const time_t&);
void thermalization_analysis(const std::string&);

#endif


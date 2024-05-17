#include <random>
#include <cmath>
#include <complex>

using std::vector;

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
vector<double> launder(vector<long int> histPoints, double min, double max, int length);

vector<long int> binCounts(vector<double> data, double min, double max, double binWidth, int length);

#endif
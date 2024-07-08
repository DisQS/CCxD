#include <random>
#include <cmath>
#include <complex>

using std::vector;

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
vector<long double> launder(vector<long double> histPoints, long double min, long double max, int length, long double binWidth, mt19937_64 gen);

vector<long double> binCounts(vector<long double> data, long double min, long double max, long double binWidth, int length);

#endif
#include <random>
#include <cmath>
#include <complex>

using std::vector;

#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
vector< double> launder(vector< double> histPoints,  double min,  double max, int length,  double binWidth,  randNums RNG);

vector< double> binCounts(vector< double> data,  double min,  double max,  double binWidth, int length);

#endif
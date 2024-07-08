#include <random>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>
// Conditional clause importing Eigen from either usr/include or from local
#if __has_include("../Eigen/eigen-master/Eigen/Dense")
#include "../Eigen/eigen-master/Eigen/Dense"
#else
#include <Eigen/Dense>
#endif


using std::sin;
using std::cos;
using std::exp;
using std::acos;


using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXcd;
const std::complex<long double> i(0.0,1.0);
const long double twopi = acos(0.0) * 4;

Matrix<std::complex<long double>,20,20> matrixReturnTRI(vector<long double> p, vector<long double> t, vector<long double> x);
Matrix<std::complex<long double>,20,1> inputVectorReturnTRI(vector<long double> p,vector<long double> t,vector<long double> inputs);
double renormalise(vector<long double> angleVector, vector<long double> scatteringAngleVector, vector<long double> phases, vector< long double> inputs);
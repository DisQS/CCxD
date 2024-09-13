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
using namespace std;

using std::sin;
using std::cos;
using std::exp;
using std::acos;


using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXcd;
const std::complex< double> i(0.0,1.0);
const  double twopi = acos(0.0) * 4;

Matrix<std::complex< double>,20,20> matrixReturnTRI(vector< double> p, vector< double> t, vector< double> x);
Matrix<std::complex< double>,20,1> inputVectorReturnTRI(vector< double> p,vector< double> t,vector< double> inputs);
Matrix<std::complex< double>,10,10> matrixReturnORIGINALT(vector< double> t, vector< double> r, vector< double> x);
Matrix<std::complex< double>,10,1> inputVectorReturnORIGINALT(vector< double> t,vector< double> r,vector< double> inputs);

double renormaliseTRI(vector< double> angleVector, vector< double> scatteringAngleVector, vector< double> phases, vector<  double> inputs);
double renormaliseORIGINALT(vector<double> tvals, vector<double> rvals, vector<double> phases, vector<double> inputs);
double renormaliseORIGINALTANALYTIC(vector<double> tvals, vector<double> rvals, vector<double> phases, vector<double> inputs);
double renormaliseORIGINALTNOINVERSE(vector<double> tvals, vector<double> rvals, vector<double> phases, vector<double> inputs);
double renormaliseORIGINALTNOINVERSEFASTERMAYBE(vector<double> tvals, vector<double> rvals, vector<double> phases, vector<double> inputs);
double renormaliseTRIINVERSE(vector<double> angleVector,  vector<double> scatteringAngleVector, vector<double> phases, vector<double> inputs);

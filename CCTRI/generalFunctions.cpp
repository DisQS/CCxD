#include "TRIRGConfig.h"
#include "randNumGen.h"
#include "histogramOperations.h"
#include "RGStep.h"
#include <iostream>

// Conditional clause importing Eigen from either usr/include or from local
#if __has_include("<Eigen/Dense>")
# include <Eigen/Dense>
#elif __has_include("../Eigen/eigen-master/Eigen/Dense")
#include "../Eigen/eigen-master/Eigen/Dense"
#endif

#include <cmath>
#include <complex>
#include <random>
#include <chrono>
#include <string>
#include <fstream>
#include <limits.h>
#include <unistd.h>
using namespace std::chrono;



using std::sin;
using std::cos;
using std::exp;
using std::acos;

 
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXcd;
//const std::complex<double> i(0.0,1.0);
//const double twopi = acos(0.0) * 4;
using std::vector;
using std::rand;
using std::mt19937_64;

const double seed = SEED;
mt19937_64 re(seed);

randNums randNums;



/*
---------------------------
randDouble (array)
---------------------------

usage:
    supplies an array of random number between the bounds provided

parameters:
    lower bound
    upper bound
    length

returns:
    an array of random doubles with predefined length

*/

vector<double> randNums::randDouble(double lower,double upper, int length) {
    std::uniform_real_distribution<double> unif(lower,upper);
    vector<double> r(length);
    for(int i{0};i<length;i++){
        r[i] = unif(re);
    }
    return r;
}
/*
---------------------------
randDouble
---------------------------

usage:
    supplies a random number between the bounds provided

parameters:
    lower bound
    upper bound

returns:
    a random double between the two given numbers

*/
double randNums::randDouble(double lower,double upper){
    std::uniform_real_distribution<double> unif(lower,upper);
    double r = unif(re);
    return r;
}

/*
---------------------------
randInt (array)
---------------------------

usage:
    generates an array of random integers between the bounds provided

parameters:
    integer lower bound
    integer upper bound
    integer length of array

returns:
    an array of random integers with predefined length
*/
vector<int> randNums::randInt(int lower, int upper, int length){
    std::uniform_int_distribution<int> unif(lower,upper);
    vector<int> r(length);
    for(int i{0};i<length;i++){
        r[i] = unif(re);
    }
    return r;
}

/*
---------------------------
randInt
---------------------------

usage:
    generates a random integer between the bounds provided

parameters:
    integer lower bound
    integer upper bound

returns:
    an random integer
*/
int randNums::randInt(int lower,int upper){
    std::uniform_int_distribution<int> unif(lower,upper);
    int r = unif(re);
    return r;
}



 
/*
---------------------------
matrixReturnTRI
---------------------------

usage:
    given initial parameters, returns the matrix representing a TRI RG unit for those parameters

parameters:
    values for node spin-mixing angles p1 through 5 (0<=px<2pi)
    normal scattering angles t  (0<=t[x]<2pi)
    8 phases x representing phase accrued over one of the 8 edges in the RG node

returns:
    20x20 (consisting of complex doubles) matrix which represents the TRI RG cell entirely

side effects:
    none

*/
Matrix<std::complex<double>,20,20> matrixReturnTRI(vector<double> p, vector<double> t, vector<double> x){
        Matrix<std::complex<double>,20,20> r {
                {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cos(p[0])*sin(t[0])*exp(-i*x[0]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sin(p[0])*exp(i*x[1]), 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0, 0.0, -i*cos(p[0]) *cos(t[0]) *exp(i*x[2]), 0.0, 0.0, sin(p[0])*exp(-i*x[0]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cos(p[0])*sin(t[0])*exp(i*x[1]), 0.0, 0.0},
                {0.0, 0.0, 1.0, 0.0, 0.0, -cos(p[0])*sin(t[0]) * exp(i*x[2]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -i * cos(p[0]) * cos(t[0]) * exp(i*x[1]), 0.0, 0.0},
                {0.0, 0.0, 0.0, 1.0, 0.0, sin(p[0]) * exp(i *x[2]), 0.0, 0.0, -i*cos(p[0]) * cos(t[0]) * exp(-i*x[0]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
 
                {-i* cos(p[1]) *cos(t[1]) *exp(-i*x[2]), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,  0.0, -sin(p[1])*exp(i *x[3]), 0.0, -cos(p[1]) *sin(t[1])* exp(-i *x[4]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -cos(p[1])* sin(t[1])* exp(i*x[3]), 0.0, sin(p[1]) *exp(-i*x[4]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {-sin(p[1]) *exp(-i*x[2]), 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,  0.0, -i *cos(p[1]) *cos(t[1])* exp(i* x[3]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {-cos(p[1]) *sin(t[1]) *exp(-i*x[2]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -i *cos(p[1])* cos(t[1]) *exp(-i *x[4]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
 
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cos(p[2])* sin(t[2])* exp(-i*x[3]), 1.0, 0.0, 0.0, 0.0,  0.0, -sin(p[2])* exp(i*x[5]), 0.0, 0.0, -i* cos(p[2])* cos(t[2])* exp(-i *x[7]), 0.0, 0.0, 0.0},
                {0.0, 0.0, -i *cos(p[2]) *cos(t[2]) *exp(i*x[0]), 0.0, 0.0, 0.0, 0.0, sin(p[2])* exp(-i*x[3]), 0.0, 1.0, 0.0, 0.0, 0.0, -cos(p[2]) *sin(t[2]) *exp(i*x[5]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, -cos(p[2])* sin(t[2]) *exp(i*x[0]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,  0.0, -i *cos(p[2]) *cos(t[2]) *exp(i*x[5]), 0.0, 0.0, -sin(p[2])* exp(-i*x[7]), 0.0, 0.0,  0.0},
                {0.0, 0.0, sin(p[2]) *exp(i*x[0]), 0.0, 0.0, 0.0,  0.0, -i *cos(p[2])* cos(t[2])* exp(-i *x[3]), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -cos(p[2])* sin(t[2])* exp(-i*x[7]), 0.0, 0.0, 0.0},
 
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -i *cos(p[3]) *cos(t[3]) *exp(-i*x[5]), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -sin(p[3])* exp(i*x[6]), 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -i*cos(p[3])* cos(t[3])* exp(i*x[4]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -cos(p[3]) *sin(t[3])* exp(i*x[6]), 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cos(p[3])* sin(t[3])* exp(i*x[4]), 0.0, 0.0, 0.0, 0.0, -sin(p[3])* exp(-i*x[5]), 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -i *cos(p[3])* cos(t[3]) *exp(i*x[6]), 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sin(p[3]) *exp(i*x[4]), 0.0, 0.0, 0.0, 0.0, -cos(p[3]) *sin(t[3]) *exp(-i*x[5]), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
 
                {0.0, 0.0, 0.0, -i *cos(p[4])* cos(t[4]) *exp(-i*x[1]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -cos(p[4]) *sin(t[4])* exp(-i*x[6]), 1.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -i *cos(p[4]) *cos(t[4]) *exp(i*x[7]), 0.0, 0.0, 0.0, 0.0, 0.0, sin(p[4])* exp(-i*x[6]), 0.0, 1.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, -sin(p[4]) *exp(-i*x[1]), 0.0, 0.0, 0.0, 0.0,  0.0, -cos(p[4]) *sin(t[4]) *exp(i*x[7]), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                {0.0, 0.0, 0.0, -cos(p[4])* sin(t[4]) *exp(-i*x[1]), 0.0, 0.0, 0.0, 0.0, 0.0, sin(p[4]) *exp(i*x[7]), 0.0, 0.0, 0.0, 0.0, 0.0, -i *cos(p[4])* cos(t[4])* exp(-i*x[6]), 0.0, 0.0, 0.0, 1.0},
        };
        return r;


}
/*
---------------------------
inputVectorReturnTRI
---------------------------

usage:
    given initial parameters, returns the vector representing the input channels of the TRI unit cell

parameters:
    values for node spin-mixing angle p (0<=p[x]<2pi)
    normal scattering angles t (0<=t[x]<2pi)
    the 4 possible inpout channels (0<=inputs[x]<1)

returns:
    1x20 (consisting of complex doubles) matrix/vector representing the inputs of the TRI unit cell
*/
Matrix<std::complex<double>,20,1> inputVectorReturnTRI(vector<double> p,vector<double> t,vector<double> in){
    Matrix<std::complex<double>,20,1> result {
        i* cos(p[0]) *cos(t[0]) *in[0],
        0,
        sin(p[0]) *in[0],
        cos(p[0])* sin(t[0]) *in[0],

        0,
        i*cos(p[1]) *cos(t[1])* in[3],
        cos(p[1])* sin(t[1])*in[3],
        -sin(p[1])* in[3],

        0, 0, 0, 0,

        cos(p[3])*sin(t[3]) *in[1],
        -sin(p[3])*in[1],
        0,
        i*cos(p[3]) *cos(t[3])*in[1],

        sin(p[4])*in[2],
        cos(p[4])* sin(t[4])* in[2],
        i * cos(p[4])* cos(t[4])*in[2],
        0
        };
    return result;
}





/*
---------------------------
launder
---------------------------

usage:
    given a histogram, it generates a distribution according to that histogram (I call it launder because it is like laundering the points in the distribution)

parameters:
    histPoints: array of histogram points
    min: minimum value that can be attained in the distribution
    max: maximum value that can be attained in the distribution
    length: how many values to generate

returns:
    a vector of specified length consisting of values that conform to the histogram provided

*/
vector<double> launder(vector<long int> histPoints, double min, double max, int length, double binWidth){
    double histPointsSum = std::reduce(histPoints.begin(),histPoints.end());
    //std::cout << histPointsSum <<std::endl;
    // make the histogram points normed by dividing each element by total 
    vector<double> normed(histPoints.size());
    for(int i{0};i<histPoints.size();i++){
        normed[i] = (histPoints[i]/histPointsSum);
    }
    // find maximum value from the histogram to create a cutoff for generating values
    double hmax  = *std::max_element(normed.begin(),normed.end());
    //std::cout << hmax << std::endl;
    // initialise laundered array
    vector<double> laundered(length);
    // populate laundered array
    //for(int i{0};i<normed.size();i++){
    //    std::cout << normed[i] << std::endl;
   // }
    for(int i{0};i<length;i++){
        // generate a 'prospect point' which could potentially be within the distribution
        vector<double> prospectPoint = {randNums.randDouble(min,max),randNums.randDouble(0,hmax)};
        // determine which bin this prospect number should fall into
        int binNo = std::floor((prospectPoint[0]-min)/binWidth);
        // keep generating prospect points until one falls within the distribution
        while(normed[binNo] < prospectPoint[1]){
            prospectPoint[0] = randNums.randDouble(min,max);
            prospectPoint[1] = randNums.randDouble(0,hmax);
            binNo = std::floor((prospectPoint[0]-min)/binWidth);

        }
        // add the point to the distribution
        laundered[i] = prospectPoint[0];

    }
    return laundered;
}

/*
---------------------------
binCounts
---------------------------

usage:
    given a set of datapoints, place them into bins of set width, returning the datapoints of the histogram representing the data

parameters:
    data: an array of data points to bin
    min: minimum value of the range of the data points
    max: maximum value of the range of the data points
    binWidth: desired width of bins

returns:
    int array

*/
vector<long int> binCounts(vector<double> data, double min, double max, double binWidth, int length){
    int amountOfBins = (int)(std::ceil((max-min)/binWidth));
    vector<long int> bins(amountOfBins);
    for(int i{0};i<length;i++){
        int binNo = (int)std::floor((data[i]-min)/binWidth);
        if(binNo >= amountOfBins){
            std::cout << binNo << std::endl;
        } else if(binNo < 0){
            std::cout << binNo <<std::endl;
        }else if(0<=binNo<amountOfBins){
            bins[binNo]+=1;
        }
    }
    return bins;
}


/*
---------------------------
renormalise
---------------------------

usage:
    performs an analytic renormalisation transformation based on an RSRG unit cell

parameters:
    angleVector: double array of spin mixing angles
    scatteringAngleVector: double array of scattering angles (thetas)
    phases: double array of random phases between nodes of the unit cell
    inputs: double array of input values into the system

returns:
    a single double representing the renormalised scattering angle for the whole unit cell

*/
double renormalise(vector<double> angleVector, vector<double> scatteringAngleVector, vector<double> phases, vector<double> inputs){
    Matrix<std::complex<double>,20,20> system = matrixReturnTRI(angleVector,scatteringAngleVector,phases);
    Matrix<std::complex<double>,20,1> inputvec = inputVectorReturnTRI(angleVector,scatteringAngleVector,inputs);
    Matrix<std::complex<double>,20,20> inv = system.inverse();
    Matrix<std::complex<double>,20,1> tmp = inv*inputvec;
    double tval = std::asin(std::abs(tmp[19])/(cos(std::asin(std::sqrt(std::pow(std::abs(tmp[14]),2) + std::pow(std::abs(tmp[1]),2))))));
    //vector<double> tval = {tmp[1],tmp[]}
    return tval;

}
#include <iostream>
#include <Eigen/Dense>
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
using std::rand;
using std::mt19937_64;
using std::vector;

 
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXcd;
const std::complex<double> i(0.0,1.0);
const double twopi = acos(0.0) * 4;
const double seed = 12;
const std::string filename = "TRIRG";
mt19937_64 re(seed);

/*

---------------------------
getPath
---------------------------

usage: 
    to create a string of the current directory the file is in

parameters:
    none
return:
    a string of the filepath
*/
std::string getPath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  std::string fullpath = std::string( result, (count > 0) ? count : 0 );
  return fullpath.substr(0,fullpath.size()-filename.size());
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
Matrix<std::complex<double>,20,1> inputVectorReturnTRI(vector<double> p,vector<double> t,vector<double> inputs){
    Matrix<std::complex<double>,20,1> result {
        i* cos(p[0]) *cos(t[0]) *inputs[0], 0, sin(p[0]) *inputs[0],cos(p[0])* sin(t[0]) *inputs[0],
         0, i*cos(p[1]) *cos(t[1])* inputs[1],cos(p[1])* sin(t[1])*inputs[1], -sin(p[1])* inputs[1],
         0, 0, 0, 0,
        cos(p[3])*sin(t[3]) *inputs[2], -sin(p[3])*inputs[2], 0, i*cos(p[3]) *cos(t[3])*inputs[2],
        sin(p[4])*inputs[3], cos(p[4])* sin(t[4])* inputs[3], i * cos(p[4])* cos(t[4])*inputs[3], 0};
    return result;
}


/*
---------------------------
randDoubleArray
---------------------------

usage:
    supplies an array of random number between the bounds provided

parameters:
    lower bound
    upper bound
    length

returns:
    an array of random doubles in an array of specified length

*/

vector<double> randDoubleArray(double lower,double upper, int length){
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
double randDouble(double lower,double upper){
    std::uniform_real_distribution<double> unif(lower,upper);
    double r = unif(re);
    return r;
}

/*
---------------------------
randIntArray
---------------------------

usage:
    generates an array of random 
*/
vector<int> randIntArray(int lower, int upper, int length){
    std::uniform_int_distribution<int> unif(lower,upper);
    vector<int> r(length);
    for(int i{0};i<length;i++){
        r[i] = unif(re);
    }
    return r;
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
vector<double> launder(vector<int> histPoints, double min, double max, int length){
    // determine how wide each bin is
    double binWidth = (max-min)/histPoints.size();
    double histPointsSum = std::reduce(histPoints.begin(),histPoints.end());
    // make the histogram points normed by dividing each element by total * binWidth
    for(double i : histPoints){
        i = i / (histPointsSum * binWidth);
    }
    // find maximum value from the histogram to create a cutoff for generating values
    double hmax  = *std::max_element(histPoints.begin(),histPoints.end());
    // initialise laundered array
    vector<double> laundered(length);
    // populate laundered array
    for(int i{0};i<length;i++){
        // generate a 'prospect point' which could potentially be within the distribution
        vector<double> prospectPoint = {randDouble(min,max),randDouble(0,hmax)};
        // determine which bin this prospect number should fall into
        int binNo = std::floor((prospectPoint[0]-min)/binWidth);
        // keep generating prospect points until one falls within the distribution
        while(histPoints[binNo] < prospectPoint[1]){
            prospectPoint = randDoubleArray(min,max,2);
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
vector<int> binCounts(double* data, double min, double max, double binWidth, int length){
    int amountOfBins = (int)std::ceil((max-min)/binWidth);
    vector<int> bins(amountOfBins);
    for(int i{0};i<length;i++){
        int binNo = (int)std::floor((data[i]-min)/binWidth);
        bins[binNo]+=1;
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
    double tval = std::asin(std::abs(tmp(19))/(cos(std::asin(std::sqrt(std::pow(std::abs(tmp(14)),2) + std::pow(std::abs(tmp(1)),2))))));
    return tval;

}

int main()
{
    // INITIALISATION PARAMS
    // -------------------------------------------
    const double angle = twopi/16; // pi/4
    vector<double> angleVector{angle,angle,angle,angle,angle};
    vector<double> inputs{1,0,0,0};
    const int length = 1000000;
    int step = 0;
    const int steps = 3;
    bool symmetrise =false;
    // -------------------------------------------

    // initial uniform distribution of theta values
    double* tprime =  (double*) malloc(length * sizeof(double));
    double* tdist = (double*) malloc(length * sizeof(double));
    double* gdist = (double*) malloc(length * sizeof(double));
    double* zdist = (double*) malloc(length * sizeof(double));
    vector<int> binsth;
    vector<int> binst;
    vector<int> binsg;
    vector<int> binsz;
    std::string path = getPath();

    for(int i{0};i<length;i++){
        tprime[i] = randDouble(0,twopi/4);
        tdist[i] = cos(tprime[i]);
        gdist[i] = std::pow(tdist[i],2);
        zdist[i] = std::log((1/gdist[i])-1);
    }
    binsth = binCounts(tprime,0,twopi/4,0.01, length);
    binst = binCounts(tdist,0,1,0.01, length);
    binsg = binCounts(gdist,0,1,0.01, length);

    std::ofstream outputth (path + "outputth" + std::to_string(0) + ".txt");
    std::ofstream outputt (path + "outputt" + std::to_string(0) + ".txt");
    std::ofstream outputg (path + "outputg" + std::to_string(0) + ".txt");
    std::ofstream outputz (path + "outputz" + std::to_string(0) + ".txt");

    //write to theta file
    for(int i{0};i<binsth.size();i++){
        outputth << binsth[i] << std::endl;
    }
    // write to t and g files
    for(int i{0};i<binst.size();i++){
        outputt << binst[i] << std::endl;
        outputg << binsg[i] << std::endl;
    }
    // write to z file
    for(int i{0};i<binsz.size();i++){
        outputz << binsz[i] << std::endl;
    }
    // close all files
    outputth.close();
    outputt.close();
    outputg.close();
    outputz.close();

    free(tdist);
    free(gdist);
    free(zdist);

    // begin clock
    auto start = high_resolution_clock::now();

    vector<int> oldTValsIndex(5);
    vector<double> oldTVals(5);

    

    
    for(int k{0};k<steps;k++){

        // initialise new array of theta values
        double* newtprime = (double*) malloc(length * sizeof(double));
        double* tdist = (double*) malloc(length * sizeof(double));
        double* gdist = (double*) malloc(length * sizeof(double));
        double* zdist = (double*) malloc(length * sizeof(double));

        // create new t values from old values
        for(int i{0};i<length;i++){
            // 5 random integers to pick the index from the old t values
            oldTValsIndex = randIntArray(0,length,5);
            for(int j{0};j<5;j++){
                oldTVals[j] = tprime[oldTValsIndex[j]];
            }
            // generate renormalised t value based on input t values, and other predefined parameters
            newtprime[i] = renormalise(angleVector,oldTVals,randDoubleArray(0,twopi,8),inputs);
            tdist[i] = cos(newtprime[i]);
            gdist[i] = std::pow(tdist[i],2);
            zdist[i] = std::log((1/gdist[i])-1);

        }
    
        // create histogram from new data

        binsz = binCounts(zdist,-25,25,0.05, length);
        if(symmetrise){
            for(int i{0};i<std::floor(binsz.size()/2);i++){
                binsz[i] = (binsz[i] + binsz[binsz.size()-i])/2;
            }
            for(int i{(int) std::floor(binsz.size()/2)};i<binsz.size();i++){
                binsz[i] = binsz[binsz.size()-i];
            }
             binsth = binCounts(newtprime,0,twopi/4,0.01, length);
            binst = binCounts(tdist,0,1,0.01, length);
            binsg = binCounts(gdist,0,1,0.01, length);
        } else{
        binsth = binCounts(newtprime,0,twopi/4,0.01, length);
        binst = binCounts(tdist,0,1,0.01, length);
        binsg = binCounts(gdist,0,1,0.01, length);
        }
        // open files to write to
        std::ofstream outputth (path + "outputth" + std::to_string(k+1) + ".txt");
        std::ofstream outputt (path + "outputt" + std::to_string(k+1) + ".txt");
        std::ofstream outputg (path + "outputg" + std::to_string(k+1) + ".txt");
        std::ofstream outputz (path + "outputz" + std::to_string(k+1) + ".txt");
        
        //write to theta file
        for(int i{0};i<binsth.size();i++){
            for(int j{0};j<std::floor(binsth[i]/100);j++){
                std::cout << "*";
            }
            outputth << binsth[i] << std::endl;
            std::cout << std::endl;
        }
        // write to t and g files
        for(int i{0};i<binst.size();i++){
            outputt << binst[i] << std::endl;
            outputg << binsg[i] << std::endl;
        }
        // write to z file
        for(int i{0};i<binsz.size();i++){
            outputz << binsz[i] << std::endl;
        }
        // close all files
        outputth.close();
        outputt.close();
        outputg.close();
        outputz.close();
        for(int i{0};i<length;i++){
            tprime[i] = newtprime[i];
        }
        free(newtprime);
        free(tdist);
        free(gdist);
        free(zdist);
    
    }
    
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    // std::cout << inv <<std::endl;
    //std::cout << system1 << std::endl;
    std::cout << duration.count() << std::endl;
    
    return 0;
}
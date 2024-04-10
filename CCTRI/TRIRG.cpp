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
const std::string filename = "TRIRG";
randNums RNG;

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


int main(int argc, char* argv[])
{
    if (argc < 3) {
        // report version
        std::cout << argv[0] << " Version " << TRIRG_VERSION_MAJOR << "."
              << TRIRG_VERSION_MINOR << std::endl;
        std::cout << "Usage: " << argv[0] << " length   steps" << std::endl;
        return 1;
    }
    // INITIALISATION PARAMS
    // -------------------------------------------
    const double angle = twopi/16; // pi/4
    vector<double> angleVector{angle,angle,angle,angle,angle};
    vector<double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int length = std::pow(10,std::stoi(argv[1]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = std::stoi(argv[2]);
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
        tprime[i] = RNG.randDouble(0,twopi/4);
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
            oldTValsIndex = RNG.randInt(0,length,5);
            for(int j{0};j<5;j++){
                oldTVals[j] = tprime[oldTValsIndex[j]];
            }
            // generate renormalised t value based on input t values, and other predefined parameters
            newtprime[i] = renormalise(angleVector,oldTVals,RNG.randDouble(0,twopi,8),inputs);
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
            for(int j{0};j<std::floor(binsth[i]);j++){
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

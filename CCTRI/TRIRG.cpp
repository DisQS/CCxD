#include "TRIRGConfig.h"
#include "randNumGen.h"
#include "histogramOperations.h"
#include "RGStep.h"
#include <iostream>
#include <filesystem>

// Conditional clause importing Eigen from either usr/include or from local
#if __has_include("../Eigen/eigen-master/Eigen/Dense")
#include "../Eigen/eigen-master/Eigen/Dense"
#else
#include <Eigen/Dense>
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
using std::sqrt;

 
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::VectorXcd;
//const std::complex<double> i(0.0,1.0);
//const double twopi = acos(0.0) * 4;
const std::string filename = "TRIRG";
namespace fs = std::filesystem;
randNums RNG;
/*
Options:
read in BOOL
read in address STR
number of samples INT
number of steps INT
symmetrise BOOL
offsetval DOUBLE
spinangle DOUBLE


*/







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
    /*argument numbers
    noOfSamples -> argv[1]
    noOfSteps -> argv[2]
    offsetVal -> argv[3]
    spinangle -> argv[4]
    symmetrise -> argv[5]
    readIn -> argv[6]
    readInAddress -> argv[7]
    */
    if (argc < 7) {
        // report version
        std::cout << argv[0] << " Version " << TRIRG_VERSION_MAJOR << "."
              << TRIRG_VERSION_MINOR << std::endl;
        std::cout << "Usage: " << argv[0] << " noOfSamples (int)   noOFSteps (int)  offsetVal (double) spinAngle (as in 2pi/spinAngle, double)    symmetrise (bool)   readIn (bool)   readInAddress (str)" << std::endl;
        return 1;
    }
    // INITIALISATION PARAMS
    // -------------------------------------------
    const double zbound = 25;
    const double angleInput = std::stod(argv[4]);
    const double angle = twopi / std::stod(argv[4]);
    vector<double> angleVector{angle,angle,angle,angle,angle};
    vector<double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int lengthInput = std::stoi(argv[1]);
    const int length = std::pow(10,std::stoi(argv[1]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = std::stoi(argv[2]);
    bool symmetrise =std::stoi(argv[5]);
    const double offsetVal = std::stod(argv[3]);
    bool readIn = std::stoi(argv[6]);
    std::string readInAddress = argv[7];
    if(readIn){
    }else{
        //const std::chrono::time_point now{std::chrono::system_clock::now()};
 
        //const std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)};
        //fs::create_directories("./" + ymd.year() + "-" + ymd.month() + "-" + ymd.day() + "_");
    }
    
    // -------------------------------------------

    // initial uniform distribution of theta values
    //double* tprime =  (double*) malloc(length * sizeof(double));
    //double* tdist = (double*) malloc(length * sizeof(double));
    //double* gdist = (double*) malloc(length * sizeof(double));
    vector<double> zdist(length);
    vector<double> tdist(length);
    vector<double> tprime(length);
    vector<double> gdist(length);
    //vector<double> zdist = malloc(length * sizeof(double));
    vector<int> binsth;
    vector<int> binst;
    vector<int> binsg;
    vector<int> binsz;
    std::string path = getPath();

    // If the user sets readIn to 1, then an ifstream is opened at the readInAddress (which should be a z distribution)
    if(readIn){
        std::ifstream currzdist;
        currzdist.open(readInAddress);
        int element;
        int i=0;
        while(currzdist >> element){
            binsz[i++] = element;
        }
        currzdist.close();
        //Launder is invoked to create a set of samples from the distribution data that was read in
        zdist = launder(binsz,-zbound,zbound,length);
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+std::exp(zdist[i]));
            tdist[i] = std::sqrt(gdist[i]);
            tprime[i] = std::acos(tdist[i]);
        }
        //Otherwise standard initialised data is produced, set here to be a uniform distribution in theta.
    }else{
                for(int i{0};i<length;i++){
            tprime[i] = RNG.randDouble(0,twopi/4);
            tdist[i] = cos(tprime[i]);
            gdist[i] = std::pow(tdist[i],2);
            zdist[i] = std::log((1/gdist[i])-1);
        }
    }
    binsth = binCounts(tprime,0,twopi/4,0.01, length);
    binst = binCounts(tdist,0,1,0.01, length);
    binsg = binCounts(gdist,0,1,0.01, length);

    //create directories to save data to

    //fs::current_path(fs::temp_directory_path());
    fs::create_directories("./CCTRI-"+ std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput));
    for(int i{0};i<steps;i++){
        fs::create_directory("./CCTRI-"+ std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput) + "/" + std::to_string(i));
    }


    // Preparing ofstreams ot read out data into relevant files
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

    //free(tdist);
    //free(gdist);
    //free(zdist);

    // begin clock for benchmarking purposes
    auto start = high_resolution_clock::now();

    vector<int> oldTValsIndex(5);
    vector<double> oldTVals(5);

    

    
    for(int k{0};k<steps;k++){

        // initialise new array of theta values
        vector<double> newtprime(length);
        vector<double> tdist(length);
        vector<double> gdist(length);
        vector<double> zdist(length);
       // double* newtprime = (double*) malloc(length * sizeof(double));
       // double* tdist = (double*) malloc(length * sizeof(double));
       // double* gdist = (double*) malloc(length * sizeof(double));
       // double* zdist = (double*) malloc(length * sizeof(double));

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
            // Each z distribution value in the first half is taken to be the arithmetic mean of that value and the corresponding value at the other end of the distribution
            for(int i{0};i<std::floor(binsz.size()/2);i++){
                binsz[i] = (binsz[i] + binsz[binsz.size()-i])/2;
            }
            // Then each +ve z distribution value is set to be the corresponding -ve distribution value, to ensure symmetry
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
        std::ofstream outputth (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput) + "/" + std::to_string(k)+ "/thdist.txt");
        std::ofstream outputt (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput) + "/" + std::to_string(k)+ "/tdist.txt");
        std::ofstream outputg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput) + "/" + std::to_string(k)+ "/gdist.txt");
        std::ofstream outputz (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string(angleInput) + "/" + std::to_string(k)+ "/zdist.txt");
        
        //write to theta file
        std::cout << k << std::endl;
        for(int i{0};i<binsth.size();i++){
            
            //Print a histogram made of stars to the command line
            /*
            for(int j{0};j<std::floor(binsth[i]);j++){
                std::cout << "*";
            }
            */
            outputth << binsth[i] << std::endl;
            //std::cout << std::endl;
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
        //free(newtprime);
        //free(tdist);
        //free(gdist);
        //free(zdist);
    
    }
    
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    // std::cout << inv <<std::endl;
    //std::cout << system1 << std::endl;
    std::cout << duration.count() << std::endl;
    
    return 0;
}

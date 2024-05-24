#include "TRIRGConfig.h"
#include "randNumGen.h"
#include "histogramOperations.h"
#include "RGStep.h"
#include <iostream>
#include <filesystem>
#include <algorithm>

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
// if DEBUG_MODE is activated (can only be activated from modifying code) terminal readouts are activated to identify current process in action.
const int DEBUG_MODE = 1;
const int WRITE_OUT_RAW=0;
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
    singlestartingthvalue -> argv[4]
    spinangle -> argv[4]
    symmetrise -> argv[5]
    readIn -> argv[6]
    readInAddress -> argv[7]
    "CCTRI-[length-input]-[steps]-[angleInput]/[finalstep]"
    */
    if (argc < 8) {
        // report version
        std::cout << argv[0] << " Version " << TRIRG_VERSION_MAJOR << "."
              << TRIRG_VERSION_MINOR << std::endl;
        std::cout << "Usage: " << argv[0] << " noOfSamples (int)   noOFSteps (int)  offsetVal (double) multiply/divide (bool) spinAngle (as in 2pi/spinAngle, double)    symmetrise (bool)   readIn (bool)   readInAddress (str)" << std::endl;
        return 1;
    }
    // INITIALISATION PARAMS
    // -------------------------------------------
    const double zbound = 25;
    const double zbinsize = 0.1;
    const double thgtbinsize = 0.01;
    const double angleInput = std::stod(argv[5]);
    const double angle = 0.01 * twopi * (std::stod(argv[5])/2);
    const double singleAngleInput = std::stod(argv[4]);
    const double singleThValue = (twopi/2) *0.01 * std::stod(argv[4]);
    vector<double> angleVector{angle,angle,angle,angle,angle};
    vector<double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int lengthInput = std::stoi(argv[1]);
    const int length = std::pow(10,std::stoi(argv[1]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = std::stoi(argv[2]);
    bool symmetrise =std::stoi(argv[6]);
    const double offsetVal = std::stod(argv[3]);
    bool readIn = std::stoi(argv[7]);
    std::string readInAddress = argv[8];
    if(readIn){
    }else{
        //const std::chrono::time_point now{std::chrono::system_clock::now()};
 
        //const std::chrono::year_month_day ymd{std::chrono::floor<std::chrono::days>(now)};
        //fs::create_directories("./" + ymd.year() + "-" + ymd.month() + "-" + ymd.day() + "_");
    }
    
    // -------------------------------------------

    // initial uniform distribution of theta values
    vector<double> zdist(length);
    vector<double> tdist(length);
    vector<double> thdist(length);
    vector<double> gdist(length);

    vector<long int> binsth;
    vector<long int> binst;
    vector<long int> binsg;
    vector<long int> binsz;

    vector<double> testdist(100000);
    vector<double> laundertest(100000);
    vector<long int> testbins(158);
    if(DEBUG_MODE){
        std::cout << "distribution and histogram arrays initialised, fetching data path" << std::endl;
    }
    std::string path = getPath();
    if(DEBUG_MODE){
        std::cout << "file path found!" <<std::endl;
    }

    // If the user sets readIn to 1, then an ifstream is opened at the readInAddress (which should be a z distribution)
    if(readIn){
        if(DEBUG_MODE){
            std::cout<< "read-in acknowledged, opening input stream" <<std::endl;
        }
        std::ifstream currzdist;
        currzdist.open(path + "/Data/" + readInAddress + +"/zdist" + ".txt");
        int element;
        int i=0;
        while(currzdist >> element){
            binsz[i++] = element;
        }
        currzdist.close();
        if(DEBUG_MODE){
            std::cout << "histogram file successfully read in, now generating new distribution based on the histogram file" << std::endl;
        }
        //Launder is invoked to create a set of samples from the distribution data that was read in
        zdist = launder(binsz,-zbound,zbound,length,zbinsize);
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+std::exp(zdist[i]));
            tdist[i] = std::sqrt(gdist[i]);
            thdist[i] = std::acos(tdist[i]);
        }
        if(DEBUG_MODE){
            std::cout << "Distibutions created successfully" <<std::endl;
        }
        //Otherwise standard initialised data is produced, set here to be a uniform distribution in theta.
    }else{
            if(DEBUG_MODE){
                std::cout << "Creating random distribution from scratch, uniform in theta" <<std::endl;
            }
            for(int i{0};i<length;i++){
                if(singleThValue != 0){
                    thdist[i] = (twopi/4) * singleThValue;
                } else{
                    thdist[i] = RNG.randDouble(0,twopi/4);
                }
                tdist[i] = cos(thdist[i]);
                gdist[i] = std::pow(tdist[i],2);
                zdist[i] = std::log((1/gdist[i])-1);
                
        }
        if(DEBUG_MODE){
                    std::cout << "Distributions successfully created!" <<std::endl;
                }
    }
    //**************************************************
    //  LAUNDER TEST
    //**************************************************
    if(DEBUG_MODE){
        std::cout << "Testing Launder function" <<std::endl;
    }
    for(int i{0};i<100000;i++){
        testdist[i] = RNG.randDouble(0,twopi/4);
        //if(i % 10 == 0){
          //  std::cout << testdist[i] << std::endl;
        //}
    }
    for(int j{0};j<10;j++){
    testbins = binCounts(testdist,0,twopi/4,0.01,100000);
    laundertest = launder(testbins,0,1,100000,0.01);
    for(int i{0};i<10000;i++){
        testdist[i] = laundertest[i];
    }
    }

    if(DEBUG_MODE){
        std::cout << "Test Passed!" <<std::endl;
    }
    //************************************************** 
    //  APPLY OFFSETVAL
    //**************************************************
    if(offsetVal != 0){
        for(int i{0};i<length;i++){
            zdist[i] = zdist[i] + offsetVal;
            gdist[i] = 1/(std::exp(zdist[i])+1);
            tdist[i] =std::sqrt(gdist[i]);
            thdist[i] = std::acos(tdist[i]);
        }
        
    }
    //vector<double> t(length);
    //for(int i{0};i<length;i++){
    //    t[i] = RNG.randDouble(0,twopi/4);
   // }
    //vector<long int> bint;
    //bint = binCounts(t,0,twopi/4,0.01,length);

    //**************************************************
    //  CREATING BINS
    //
    if(DEBUG_MODE){
        std::cout << "Counting bins and creating histograms from created distributions.." <<std::endl;
    }
    binsth = binCounts(thdist,0,twopi/4,thgtbinsize, length);
    binst = binCounts(tdist,0,1,thgtbinsize, length);
    binsg = binCounts(gdist,0,1,thgtbinsize, length);
    binsz = binCounts(zdist,-zbound,zbound,zbinsize,length);
    if(DEBUG_MODE){
        std::cout << "..Done!" <<std::endl <<std::endl;
    }
    //create directories to save data to

    //fs::current_path(fs::temp_directory_path());
    if(DEBUG_MODE){
        std::cout << "Creating directories.." <<std::endl;
    }
    fs::create_directories("./Data/CCTRI-"+ std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput)+std::to_string((int)singleAngleInput));
    for(int i{0};i<steps+1;i++){
        fs::create_directory("./Data/CCTRI-"+ std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(i));
    }
    if(DEBUG_MODE){
        std::cout << "..Done!" <<std::endl;
    }


    // Preparing ofstreams ot read out data into relevant files
    std::ofstream outputth (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)(angleInput) +std::to_string((int)singleAngleInput) +  "/" + std::to_string(0)+"/thdist"+ ".txt");
    std::ofstream outputt (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(0)+"/tdist" +  ".txt");
    std::ofstream outputg (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(0)+"/gdist" +  ".txt");
    std::ofstream outputz (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(0)+"/zdist" + ".txt");

    //std::ofstream rawth (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawth" + std::to_string(0) + ".txt");
    //std::ofstream rawt (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawt" + std::to_string(0) + ".txt");
    //std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+"/rawg" + std::to_string(0) + ".txt");
    std::ofstream rawz (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(0)+"/rawz" +".txt");
    if(DEBUG_MODE){
        std::cout << "Writing to files.." <<std::endl;
    }
    //write to theta file
    for(int i{0};i<binsth.size();i++){
        outputth << binsth[i] << std::endl;
    //    std::cout << bint[i] << std::endl;
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
        //rawth << thdist[i] << std::endl;
        //rawt << tdist[i] << std::endl;
        //rawg << gdist[i] << std::endl;
        //rawz << zdist[i] << std::endl;
    }

    //rawth.close();
    //rawt.close();
    //rawg.close();
    rawz.close();
    if(DEBUG_MODE){
        std::cout << "..Done!" <<std::endl;
    }
    // begin clock for benchmarking purposes


    

    
    for(int k{0};k<steps;k++){
        
        auto start = high_resolution_clock::now();
        std::cout << "Starting " << std::to_string(k+1) <<"th iteration" <<std::endl;
        // initialise new array of theta values
        vector<double> oldthdist(length);
        for(int i{0};i<length;i++){
            oldthdist[i] = thdist[i];
        }
        std::cout << "Renormalising" <<std::endl;
        // create new t values from old values
        for(int i{0};i<length;i++){

            vector<int> oldTValsIndex(5);
            vector<double> oldTVals(5);
            // 5 random integers to pick the index from the old t values
            //oldTValsIndex = RNG.randInt(0,(length-1),5);
            for(int j{0};j<5;j++){
                
                oldTVals[j] = oldthdist[RNG.randInt(0,length-1)];
            }
            // generate renormalised t value based on input t values, and other predefined parameters
            thdist[i] = renormalise(angleVector,oldTVals,RNG.randDouble(0,twopi,8),inputs);
            //std::cout << thdist[i] << std::endl;
            tdist[i] = cos(thdist[i]);
            gdist[i] = std::pow(tdist[i],2);
            zdist[i] = std::log((1/gdist[i])-1);

        }
        std::cout << "Renormalised! Now creating new distributions from data" << std::endl;
        // create histogram from new data

        int zbinlength = 2 * zbound / zbinsize;
        vector<long int> newbinz(zbinlength);
        vector<long int> newbinzrev(zbinlength);
        vector<long int> newbinth;
        vector<long int> newbint;
        vector<long int> newbing;
        newbinz = binCounts(zdist,-zbound,zbound,zbinsize, length);

        if(symmetrise){
            std::cout << "The z distribution is going to be symmetrised" <<std::endl;
            // Each z distribution value in the first half is taken to be the arithmetic mean of that value and the corresponding value at the other end of the distribution
            std::reverse_copy(std::begin(newbinz),std::end(newbinz),std::begin(newbinzrev));
            vector<double> symdist(length);
            for(int i{0};i< zbinlength;i++){
                newbinz[i] = (newbinz[i] +newbinzrev[i])/2;
            }
            std::cout << "Laundering symmetrised distribution.." <<std::endl;
            symdist = launder(newbinz,-zbound,zbound,length,zbinsize);
            std::cout << "..Done!" <<std::endl;
            //for(int i{0};i<symdist.size();i++){
            //    std::cout << symdist[i] <<std::endl;
           // }
            for(int i{0};i<length;i++){
                zdist[i] = symdist[i];
                gdist[i] = 1/(std::exp(zdist[i])+1);
                tdist[i] =std::sqrt(gdist[i]);
                thdist[i] = std::acos(tdist[i]);
            }
            newbinth = binCounts(thdist,0,twopi/4,thgtbinsize, length);

            newbint = binCounts(tdist,0,1,thgtbinsize, length);

            newbing = binCounts(gdist,0,1,thgtbinsize, length);


        } else{

        newbinth = binCounts(thdist,0,twopi/4,thgtbinsize, length);
        newbint = binCounts(tdist,0,1,thgtbinsize, length);
        newbing = binCounts(gdist,0,1,thgtbinsize, length);
        }

        // open files to write to
        std::ofstream outputth (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) +std::to_string((int)singleAngleInput) +  "/" + std::to_string(k+1)+ "/thdist.txt");
        std::ofstream outputt (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(k+1)+ "/tdist.txt");
        std::ofstream outputg (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) +std::to_string((int)singleAngleInput) +  "/" + std::to_string(k+1)+ "/gdist.txt");
        std::ofstream outputz (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(k+1)+ "/zdist.txt");

       // std::ofstream rawth (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/thraw.txt");
       // std::ofstream rawt (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/traw.txt");
       // std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/graw.txt");
        std::ofstream rawz (path + "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + std::to_string((int)singleAngleInput) + "/" + std::to_string(k+1)+ "/zraw.txt");
        
        //write to theta file
        if(WRITE_OUT_RAW==1){
            for(int i{0};i<length;i++){
          //  rawth << thdist[i] << std::endl;
          //  rawt << tdist[i] << std::endl;
         //   rawg << gdist[i] << std::endl;
                //rawz << zdist[i] << std::endl;

            }
        }
        rawz.close();
        for(int i{0};i<binsth.size();i++){
            
            //Print a histogram made of stars to the command line
            /*
            for(int j{0};j<std::floor(binsth[i]);j++){
                std::cout << "*";
            }
            */
            binsth[i] = newbinth[i];
            outputth << binsth[i] << std::endl;
            //std::cout << std::endl;
        }
        // write to t and g files
        for(int i{0};i<binst.size();i++){
            binst[i] = newbint[i];
            binsg[i] = newbing[i];
            outputt << binst[i] << std::endl;
            outputg << binsg[i] << std::endl;
        }
        // write to z file
        for(int i{0};i<binsz.size();i++){
            binsz[i] = newbinz[i];
            outputz << newbinz[i] << std::endl;
        }
        // close all files
        outputth.close();
        outputt.close();
        outputg.close();
        outputz.close();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        // std::cout << inv <<std::endl;
        //std::cout << system1 << std::endl;
        std::cout << duration.count()/1000000 << std::endl;
    }

    
    
    
    
    return 0;
}

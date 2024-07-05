#include "TRIRGConfig.h"
#include "randNumGen.h"
#include "histogramOperations.h"
#include "RGStep.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <sstream>
#include <vector>
#include <iterator>
#include <mpi.h>

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
using std::mt19937_64;
//const std::complex<double> i(0.0,1.0);
//const double twopi = acos(0.0) * 4;
const std::string filename = "TRIRG";
namespace fs = std::filesystem;
// if DEBUG_MODE is activated (can only be activated from modifying code) terminal readouts are activated to identify current process in action.


double seed = SEED;


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
//*********************************
// GLOBAL OPTIONS
//*********************************
// DEBUG_MODE outputs more to console if set to 1
const int DEBUG_MODE = 1;
// SEP_FILE_OUTPUT outputs distributions into separate files if set to 1
const int SEP_FILE_OUTPUT = 0;
// OUTPUT_FREQ outputs every nth distribution when set to n
const int OUTPUT_FREQ = 4;
// WRITE_OUT_RAW writes out the raw z data if set to 1
const int WRITE_OUT_RAW=0;
//**********************************






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


template <typename Out>
void split(const std::string &s, char delim, Out result) {
    std::istringstream iss(s);
    std::string item;
    while (std::getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}


int main(int argc, char* argv[])
{

    int current_rank;
    int initialisation_error;
    int master;
    int number_of_processes;
    MPI_Status status;
    //MPI INITIALISATION

    initialisation_error = MPI_Init(&argc,&argv);
    if(initialisation_error != 0){
        std::cout << "\n";
        std::cout << "Fatal error!\n";
        std::cout << "MPI_Init returned ierr = " << initialisation_error << "\n";
        exit(1);
    }

    //Get the number of processes

    MPI_Comm_size(MPI_COMM_WORLD,&number_of_processes);

    //Determine rank of current process

    MPI_Comm_rank(MPI_COMM_WORLD,&current_rank);

    //std::cout << "I am rank" << current_rank <<"\n";
    
    //Each process creates their own number generator with varying seed
    mt19937_64 re(seed + current_rank);
    randNums RNG;
    //std::cout << RNG.randDouble(0.0,1.0,re) << "\n";

    /*argument numbers
    noOfSamples -> argv[1]
    noOfSteps -> argv[2]
    offsetVal -> argv[3]
    singlestartingthvalue -> argv[4]
    spinangle -> argv[5]
    symmetrise -> argv[6]
    readIn -> argv[7]
    readInAddress -> argv[8]
    "CCTRI-[length-input]-[steps]-[angleInput]/[finalstep]"
    */
    if (argc < 1) {
        // report version
        std::cout << argv[0] << " Version " << TRIRG_VERSION_MAJOR << "."
              << TRIRG_VERSION_MINOR << std::endl;
        std::cout << "Usage: " << argv[0] << " noOfSamples (int)   noOFSteps (int)  offsetVal (double) multiply/divide (bool) spinAngle (as in 2pi/spinAngle, double)    symmetrise (bool)   readIn (bool)   readInAddress (str)" << std::endl;
        std::cout << argv[0] << argv[2] << argv[2] << argv[3] << argv[4] << argv[5] << argv[6] << std::endl;
        return 1;
    }
    // INITIALISATION PARAMS
    // -------------------------------------------
    std::vector<std::string> arguments = split(argv[1],' ');
    const double zbound = 25;
    const double zbinsize = 0.01;
    const double thgtbinsize = 0.001;
    const double angleInput = std::stod(arguments[4]);
    const double angle = 0.01 * twopi * (std::stod(arguments[4])/2);
    const double singleAngleInput = std::stod(arguments[3]);
    const double singleThValue = (twopi/2) *0.01 * std::stod(arguments[3]);
    vector<double> angleVector{angle,angle,angle,angle,angle};
    vector<double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int lengthInput = std::stoi(arguments[0]);
    const int length = std::pow(10,std::stoi(arguments[0]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = std::stoi(arguments[1]);
    bool symmetrise =std::stoi(arguments[5]);
    const double offsetVal = std::stod(arguments[2]);
    int readIn = std::stoi(arguments[6]);
    std::string readInAddress = arguments[7];

    std::cout << "Z bound: " << zbound << std::endl;
    std::cout << "Z bin size: " << zbinsize << std::endl;
    std::cout << "Unit bin size: " << thgtbinsize << std::endl;
    std::cout << "Amount of samples: " << length << std::endl;
    std::cout << "Amount of steps: "  << steps << std::endl;
    std::cout << "Symmetrise set?: " << symmetrise << std::endl;
    std::cout << "Offset value: " << offsetVal << std::endl;

    if(readIn==1){
    }else{
        
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
    if(readIn==1){
        if(DEBUG_MODE){
            std::cout<< "read-in acknowledged, opening input stream" <<std::endl;
        }
        std::ifstream currzdist;
        currzdist.open(path + "/Data/" + readInAddress + +"/distsfinal" + ".txt");
        std::string element;
        int i=0;
        while(currzdist >> element){
            binsz[i++] = std::stoi(element);
        }
        currzdist.close();
        if(DEBUG_MODE){
            std::cout << "histogram file successfully read in, now generating new distribution based on the histogram file" << std::endl;
        }
        //Launder is invoked to create a set of samples from the distribution data that was read in
        zdist = launder(binsz,-zbound,zbound,length,zbinsize,re);
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
    }else if(readIn==2){
        if(DEBUG_MODE){
            std::cout<< "read-in acknowledged, opening input stream for raw data" <<std::endl;
        }
        std::ifstream currzraw;
        currzraw.open(path + "/Data/" + readInAddress + +"/rawfinal" + ".txt");
        std::string element;
        int i=0;
        while(currzraw >> element){
            zdist[i++] = std::stod(element);
        }
        currzraw.close();
        if(DEBUG_MODE){
            std::cout << "raw file successfully read in, now converting to other distributions based on the raw file" << std::endl;
        }
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+std::exp(zdist[i]));
            tdist[i] = std::sqrt(gdist[i]);
            thdist[i] = std::acos(tdist[i]);
        }

    }else{
            if(DEBUG_MODE){
                std::cout << "Creating random distribution from scratch, uniform in theta" <<std::endl;
            }
            for(int i{0};i<length;i++){
                if(singleThValue != 0){
                    thdist[i] = singleThValue;
                } else{
                    thdist[i] = RNG.randDouble(0,twopi/4,re);
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
        testdist[i] = RNG.randDouble(0,twopi/4,re);
        //if(i % 10 == 0){
          //  std::cout << testdist[i] << std::endl;
        //}
    }
    for(int j{0};j<10;j++){
    testbins = binCounts(testdist,0,twopi/4,0.01,100000);
    laundertest = launder(testbins,0,1,100000,0.01,re);
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
    //**************************************************
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
    std::string outputPath =  "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "-" + std::to_string((int)singleAngleInput) +  "/";
    //fs::current_path(fs::temp_directory_path());
    if(DEBUG_MODE){
        std::cout << "Creating directories.." <<std::endl;
    }
    fs::create_directories("." + outputPath);
    //for(int i{0};i<steps+1;i++){
    //    fs::create_directory("." + outputPath + "/" + std::to_string(i));
    //}
    if(DEBUG_MODE){
        std::cout << "..Done!" <<std::endl;
    }
    vector<long int> avgbinsth(binsth.size());
    vector<long int> avgbinst(binst.size());
    vector<long int> avgbinsg(binsg.size());
    vector<long int> avgbinsz(binsz.size());

    

    if(current_rank!= 0){
        //MPI_Send();
    } else {

    }



    std::cout << "Saving to: " + outputPath << std::endl;
    // Preparing ofstreams ot read out data into relevant files
    if(SEP_FILE_OUTPUT){
        std::ofstream outputth (path + outputPath +"thdist0"+ ".txt");
        std::ofstream outputt (path + outputPath + "tdist0" +  ".txt");
        std::ofstream outputg (path + outputPath + "gdist0" +  ".txt");
        std::ofstream outputz (path + outputPath + "zdist0" + ".txt");
    
        //std::ofstream rawth (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawth" + std::to_string(0) + ".txt");
        //std::ofstream rawt (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawt" + std::to_string(0) + ".txt");
        //std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+"/rawg" + std::to_string(0) + ".txt");
        std::ofstream rawz (outputPath +"rawz0" +".txt");
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
        if(WRITE_OUT_RAW==1){
            for(int i{0};i<length;i++){
            rawz << zdist[i] << std::endl;
            }
        }


        rawz.close();
        if(DEBUG_MODE){
            std::cout << "..Done!" <<std::endl;
        }
    }else{
        
        std::ofstream allbins (path + outputPath + "dists0.txt");
        allbins << binsth.size() << std::endl;
        allbins << binst.size() << std::endl;
        allbins << binsg.size() <<std::endl;
        allbins << binsz.size() << std::endl;
        for(int i{0};i<binsth.size();i++){
            allbins << binsth[i] << std::endl;
        }
        for(int i{0};i<binst.size();i++){
            allbins << binst[i] << std::endl;
        }
        for(int i{0};i<binsg.size();i++){
            allbins << binsg[i] << std::endl;
        }
        for(int i{0};i<binsz.size();i++){
            allbins << binsz[i] << std::endl;
        }
        allbins.close();
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
        //vector<double> oldTVals(%* length);
        //oldTVals = launder(binsth,0,twopi/4,5 * length,thgtbinsize);
        for(int i{0};i<length;i++){

            vector<int> oldTValsIndex(5);
            vector<double> oldTVals(5);
            // 5 random integers to pick the index from the old t values
            oldTValsIndex = RNG.randInt(0,(length-1),re,5);
            for(int j{0};j<5;j++){
                
                oldTVals[j] = oldthdist[RNG.randInt(0,length-1,re)];
            }
            // generate renormalised t value based on input t values, and other predefined parameters
            
            //thdist[i] = renormalise(angleVector,{oldTVals[(5* i)],oldTVals[(5* i)+1],oldTVals[(5* i)+2],oldTVals[(5* i)+3],oldTVals[(5* i)]+4},RNG.randDouble(0,twopi,8),inputs);
            thdist[i] = renormalise(angleVector,oldTVals,RNG.randDouble(0,twopi,re,8),inputs);
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
            symdist = launder(newbinz,-zbound,zbound,length,zbinsize,re);
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
        if(((k+1) % OUTPUT_FREQ == 0) || k==0){
            if(SEP_FILE_OUTPUT){
                std::cout << "Saving to: " + outputPath << std::endl;
                std::ofstream outputth (path + outputPath +  "thdist" + std::to_string(k+1) + ".txt");
                std::ofstream outputt (path + outputPath +  "tdist" + std::to_string(k+1) + ".txt");
                std::ofstream outputg (path + outputPath +  "gdist" + std::to_string(k+1) + ".txt");
                std::ofstream outputz (path + outputPath +  "zdist" + std::to_string(k+1) + ".txt");

       // std::ofstream rawth (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/thraw.txt");
       // std::ofstream rawt (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/traw.txt");
       // std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/graw.txt");
                std::ofstream rawz (path + outputPath + "zraw" + std::to_string(k+1) + ".txt");
        
        //write to theta file
                if(WRITE_OUT_RAW==1){
                    for(int i{0};i<length;i++){
          //  rawth << thdist[i] << std::endl;
          //  rawt << tdist[i] << std::endl;
         //   rawg << gdist[i] << std::endl;
                        rawz << zdist[i] << std::endl;

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
            }else{
                std::ofstream allbins (path + outputPath + "dists" + std::to_string(k+1) + ".txt");
                allbins << binsth.size() << std::endl;
                allbins << binst.size() << std::endl;
                allbins << binsg.size() <<std::endl;
                allbins << binsz.size() << std::endl;
                for(int i{0};i<binsth.size();i++){
                    binsth[i] = newbinth[i];
                    allbins << binsth[i] << std::endl;
                }
                for(int i{0};i<binst.size();i++){
                    binst[i] = newbint[i];
                    allbins << binst[i] << std::endl;
                }
                for(int i{0};i<binsg.size();i++){
                    binsg[i] = newbing[i];
                    allbins << binsg[i] << std::endl;
                }
                for(int i{0};i<binsz.size();i++){
                    binsz[i] = newbinz[i];
                    allbins << binsz[i] << std::endl;
                }
                allbins.close();
            }
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        // std::cout << inv <<std::endl;
        //std::cout << system1 << std::endl;
        std::cout << duration.count()/1000000 << std::endl;
    }

    
    
    
    
    return 0;
}

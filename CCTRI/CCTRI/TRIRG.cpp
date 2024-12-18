#include "TRIRGConfig.hpp"
#include "randNumGen.hpp"
#include "histogramOperations.hpp"
#include "RGStep.hpp"
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

using namespace std;



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
//*********************************
// GLOBAL OPTIONS
//*********************************
// DEBUG_MODE outputs more to console if set to 1
const int DEBUG_MODE = 1;
// SEP_FILE_OUTPUT outputs distributions into separate files if set to 1
const int SEP_FILE_OUTPUT = 0;
// OUTPUT_FREQ outputs every nth distribution when set to n
const int OUTPUT_FREQ = 1;
// WRITE_OUT_RAW writes out the raw z data if set to 1
const int WRITE_OUT_RAW=0;
// USE_THETA uses the theta parametrisation if it is set to 1, goes back to CC2D if set to 0
const int USE_THETA = 1;

const int INITIAL_DISTRIBUTION=1;
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

    cout << "Seed is: " << seed << endl;
    int current_rank;
    int initialisation_error;
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
    cout << number_of_processes << " <- number of processes" << endl;

    //Determine rank of current process

    MPI_Comm_rank(MPI_COMM_WORLD,&current_rank);

    std::cout << "Initialised process number: " << current_rank << std::endl;
    
    //Each process creates their own number generator with varying seed
    mt19937_64 re(seed + current_rank);
    RNG.gen = re;

    /*argument numbers
    noOfSamples -> arguments[0]
    noOfSteps -> arguments[1]
    offsetVal -> arguments[2]
    singlestartingthvalue -> arguments[4]
    spinangle -> argument[3]
    symmetrise -> arguments[5]
    readIn -> arguments[6]
    readInAddress -> arguments[7]
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
    const  double zbound = 25;
    const  double zbinsize = 0.001;
    const  double thgtbinsize = 0.001;
    const double angleInput = std::stoi(arguments[3].substr(2,5));
    //const  double angleInput = std::stod(arguments[3]);
//    const  double angle = 0.01 * twopi * (std::stod(arguments[3])/2);
    const double angle = std::stod(arguments[3]);
    const  double singleAngleInput = std::stoi(arguments[4].substr(2,5));
    //const  double singleThValue = (twopi/2) *0.01 * std::stod(arguments[4]);
    const double singleThValue = std::stod(arguments[4]);
    //const double singleThValue = 0.84743;
    vector< double> angleVector{angle,angle,angle,angle,angle};
    vector< double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int lengthInput = std::stoi(arguments[0]);
    const int length = std::pow(10,std::stoi(arguments[0]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = std::stoi(arguments[1]);
    bool symmetrise =std::stoi(arguments[5]);
    const  double offsetVal = std::stod(arguments[2]);
    int readIn = std::stoi(arguments[6]);
    std::string readInAddress = arguments[7];
    if(DEBUG_MODE && current_rank == 0){
        std::cout << "Z bound: " << zbound << std::endl;
        std::cout << "Z bin size: " << zbinsize << std::endl;
        std::cout << "Unit bin size: " << thgtbinsize << std::endl;
        std::cout << "Amount of samples: " << length << std::endl;
        std::cout << "Amount of steps: "  << steps << std::endl;
        std::cout << "Symmetrise set?: " << symmetrise << std::endl;
        std::cout << "Offset value: " << offsetVal << std::endl;
        std::cout << "Spin mixing angle: " << angle << std::endl;
        
    }

    if(readIn==1){
    }else{
        
    }
    
    // -------------------------------------------

    // initial uniform distribution of theta values
    vector< double> zdist(length);
    vector< double> tdist(length);
    vector< double> thdist(length);
    vector< double> gdist(length);

    vector< double> binsth(1/thgtbinsize);
    vector< double> binst(1/thgtbinsize);
    vector< double> binsg(1/thgtbinsize);
    vector< double> binsz((2 * zbound)/zbinsize);

    vector< double> testdist(100000);
    vector< double> laundertest(100000);
    vector< double> testbins(158);
    if(DEBUG_MODE && current_rank == 0){
        std::cout << "distribution and histogram arrays initialised, fetching data path" << std::endl;
    }
    std::string path = getPath();
    if(DEBUG_MODE && current_rank == 0){
        std::cout << "file path found!" <<std::endl;
    }

    // If the user sets readIn to 1, then an ifstream is opened at the readInAddress (which should be a z distribution)
    if(readIn==1){
        if(DEBUG_MODE && current_rank == 0){
            std::cout<< "read-in acknowledged, opening input stream" <<std::endl;
        }
        std::ifstream currzdist;
        std::cout << "Opening from" << path << "/Data/" << readInAddress << std::endl;
        currzdist.open(path + "/Data/" + readInAddress);
        std::string element;
        int z=0;
        int skips = 1;
        std::cout << "test" << std::endl;
        for(int i{0};i<3;i++){
          std::getline(currzdist,element);
          std::cout << element << std::endl;
          skips+= std::stoi(element);

        }
        for(int i{0};i<skips;i++){
          std::getline(currzdist,element);
        }

        while(currzdist >> element){
            binsz[z++] = std::stod(element);
        }
        currzdist.close();
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "histogram file successfully read in, now generating new distribution based on the histogram file" << std::endl;
        }
        for(int i{25000};i<25100;i++){
            cout << binsz[i] << endl;
        }
        //Launder is invoked to create a set of samples from the distribution data that was read in
        zdist = launder(binsz,-zbound,zbound,length,zbinsize, RNG);
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+std::exp(zdist[i]));
            tdist[i] = std::sqrt(gdist[i]);
            thdist[i] = std::acos(tdist[i]);
        }
        for(int i{0};i<100;i++){
            cout << zdist[i] << endl;
        }
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "Distibutions created successfully" <<std::endl;
        }
        //Otherwise standard initialised data is produced, set here to be a uniform distribution in theta.
    }else if(readIn==2){
        if(DEBUG_MODE && current_rank == 0){
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
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "raw file successfully read in, now converting to other distributions based on the raw file" << std::endl;
        }
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+std::exp(zdist[i]));
            tdist[i] = std::sqrt(gdist[i]);
            thdist[i] = std::acos(tdist[i]);
        }

    }else{
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "Creating random distribution from scratch, uniform in theta" <<std::endl;
        }
        for(int i{0};i<length;i++){
            if(singleThValue != 0){
                thdist[i] = singleThValue;
                tdist[i] = cos(thdist[i]);
                gdist[i] = pow(tdist[i],2);
                zdist[i] = log((1/gdist[i])-1);
            } else{
                if(INITIAL_DISTRIBUTION==0){
                    thdist[i] = RNG.randDouble(0,twopi/4);
                    tdist[i] = std::cos(thdist[i]);
                    gdist[i] = std::pow(tdist[i],2);
                    zdist[i] = std::log((1/gdist[i])-1);
                } else{
                    tdist[i] = std::sqrt(RNG.randDouble(0,1));
                    thdist[i] = std::acos(tdist[i]);
                    gdist[i] = std::pow(tdist[i],2);
                    zdist[i] = log((1/gdist[i])-1);
                }
            }
            
        }
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "Distributions successfully created!" <<std::endl;
        }
    }

    //**************************************************
    //  LAUNDER TEST
    //**************************************************
    if(DEBUG_MODE && current_rank==0){
        std::cout << "Testing Launder function" <<std::endl;
    }
    for(int i{0};i<100000;i++){
        testdist[i] = RNG.randDouble(0,twopi/4);
        //if(i % 10 == 0){
          //  std::cout << testdist[i] << std::endl;
        //}
    }
    for(int j{0};j<10;j++){
    testbins = binCounts(testdist,0.0,twopi/4,0.01,100000);
    laundertest = launder(testbins,0.0,1.0,100000,0.01, RNG);
    for(int i{0};i<10000;i++){
        testdist[i] = laundertest[i];
    }
    }

    if(DEBUG_MODE && current_rank==0){
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
       for(int i{0};i<100;i++){
            cout << zdist[i] << endl;
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
    if(DEBUG_MODE && current_rank==0){
        std::cout << "Counting bins and creating histograms from created distributions.." <<std::endl;
    }
    binsth = binCounts(thdist,0,twopi/4,thgtbinsize, length);
    binst = binCounts(tdist,0,1,thgtbinsize, length);
    binsg = binCounts(gdist,0,1,thgtbinsize, length);
    binsz = binCounts(zdist,-zbound,zbound,zbinsize,length);
    for(int i{25000};i<25100;i++){
        cout << binsz[i] << endl;
    }
    if(DEBUG_MODE && current_rank==0){
        std::cout << "..Done!" <<std::endl <<std::endl;
    }
    //create directories to save data to
    cout << to_string(offsetVal).substr(to_string(offsetVal).find(".")+1,3) << endl;
    std::string outputPath =  "/Data/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "-" + std::to_string((int)singleAngleInput) +  "/" + to_string(offsetVal).substr(to_string(offsetVal).find(".")+1,3) + "/";
    //fs::current_path(fs::temp_directory_path());
    if(DEBUG_MODE && current_rank==0){
        std::cout << "Creating directories.." <<std::endl;
    }
    if(current_rank==0){
        fs::create_directories("." + outputPath);
    }
    //for(int i{0};i<steps+1;i++){
    //    fs::create_directory("." + outputPath + "/" + std::to_string(i));
    //}
    if(DEBUG_MODE && current_rank==0){
        std::cout << "..Done!" <<std::endl;
    }
    vector< double> allbins(binsth.size()+binst.size()+binsg.size()+binsz.size());
    vector< double> avgbins(binsth.size()+binst.size()+binsg.size()+binsz.size());

    //pack up bins into one array to send
    for(int i{0};i<binsth.size();i++){
        allbins[i] = binsth[i];
    }
    for(int i{0};i<binst.size();i++){
        allbins[i+binsth.size()] = binst[i];
    }
    for(int i{0};i<binsg.size();i++){
        allbins[i+binsth.size()+binst.size()] = binsg[i];
    }
    for(int i{0};i<binsz.size();i++){
        allbins[i+binsth.size()+binst.size()+binsg.size()] = binsz[i];
    }

    
    //send bins from each process to the master process
    if(current_rank!= 0){
        std::cout << "Sending bins from rank: " << current_rank << std::endl;
        MPI_Send(&allbins[0],allbins.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
    } else {
        for(int i{0};i<allbins.size();i++){
            avgbins[i] = allbins[i];
        }
        for(int i{1};i<number_of_processes;i++){
            std::cout << "Receiving bins to rank: "<< current_rank << std::endl;
            MPI_Recv(&allbins[0],allbins.size(),MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
            for(int j{0};j<allbins.size();j++){
                avgbins[j]+=allbins[j];
            }
        }
        for(int j{0};j<allbins.size();j++){
            allbins[j] = avgbins[j]/number_of_processes;
        }
    }
    //Write out and then broadcast
    if(current_rank==0){
        std::cout << "Saving to: " + outputPath << std::endl;
        // Preparing ofstreams ot read out data into relevant files
        if(SEP_FILE_OUTPUT){
            for(int i{0};i<binsth.size();i++){
                binsth[i] = allbins[i]; 
            }
            for(int i{0};i<binst.size();i++){
                binst[i] = allbins[i+binsth.size()]; 
            }
            for(int i{0};i<binsg.size();i++){
                binsg[i] = allbins[i+binsth.size()+binst.size()]; 
            }
            for(int i{0};i<binsth.size();i++){
                binsz[i] = allbins[i+binsth.size()+binst.size()+binsg.size()]; 
            }
            std::ofstream outputth (path + outputPath +"thdist0"+ ".txt");
            std::ofstream outputt (path + outputPath + "tdist0" +  ".txt");
            std::ofstream outputg (path + outputPath + "gdist0" +  ".txt");
            std::ofstream outputz (path + outputPath + "zdist0" + ".txt");
    
            //std::ofstream rawth (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawth" + std::to_string(0) + ".txt");
            //std::ofstream rawt (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawt" + std::to_string(0) + ".txt");
            //std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+"/rawg" + std::to_string(0) + ".txt");
            std::ofstream rawz (outputPath +"rawz0" +".txt");
            if(DEBUG_MODE && current_rank==0){
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
            if(DEBUG_MODE && current_rank==0){
                std::cout << "..Done!" <<std::endl;
            }
        }else{
        
            std::ofstream allbinsout (path + outputPath + "dists0.txt");
            allbinsout << binsth.size() << std::endl;
            allbinsout << binst.size() << std::endl;
            allbinsout << binsg.size() <<std::endl;
            allbinsout << binsz.size() << std::endl;
            for(int i{0};i<allbins.size();i++){
                allbinsout << allbins[i] << std::endl;
            }
            allbinsout.close();
        }

    } 

    

    
    
    



    
    // begin clock for benchmarking purposes
    

    

    
    for(int k{0};k<steps;k++){
        MPI_Bcast(&allbins[0],allbins.size(),MPI_DOUBLE,0,MPI_COMM_WORLD);
        //std::cout << std::setprecision(80) <<  allbins[0] << " " << current_rank << std::endl;
        //std::cout << allbins.size() << " " << current_rank<<std::endl;
        //std::cout << binsth.size() << " " << current_rank<<std::endl;
        //std::cout << std::setprecision(80)<<  binsth[0] << " " << current_rank<< std::endl;
        //binsth[0] = allbins[0];
        //unpack into separate bins
        for(int i{0};i<binsth.size();i++){
            binsth[i] = allbins[i]; 
        }
        for(int i{0};i<binst.size();i++){
            binst[i] = allbins[i+binsth.size()]; 
        }
        for(int i{0};i<binsg.size();i++){
            binsg[i] = allbins[i+binsth.size()+binst.size()]; 
        }
        for(int i{0};i<binsth.size();i++){
            binsz[i] = allbins[i+binsth.size()+binst.size()+binsg.size()]; 
        }
        for(int i{25000};i<25100;i++){
            cout << binsz[i] << endl;
        }
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "Starting " << std::to_string(k+1) <<"th iteration" <<std::endl;
        }
        // initialise new array of theta values
        vector< double> oldzdist(length);
        vector< double> oldthdist(length);
        oldzdist = launder(binsz,-zbound,zbound,length,zbinsize, RNG);
        if(USE_THETA==1){
            for(int i{0};i<length;i++){
                oldthdist[i] = acos(sqrt(1/(exp(oldzdist[i])+1)));
            }
        }else{    
            for(int i{0};i<length;i++){
                oldthdist[i] = (std::sqrt(1/(std::exp(oldzdist[i])+1)));
            }
        }
        if(DEBUG_MODE && current_rank == 0){
            std::cout << "Renormalising" <<std::endl;
        }
        // create new t values from old values
        //vector<double> oldTVals(%* length);
        //oldTVals = launder(binsth,0,twopi/4,5 * length,thgtbinsize);
        //vector<int> oldTValsIndex(5);
        vector< double> oldTVals(5);
        vector<double> oldRVals(5);
        vector<double> tvalstest(5);
        vector<double> rvalstest(5);
        //vector<high_resolution_clock::time_point> starts(length);
        //vector<high_resolution_clock::time_point> stops(length);
        

        auto start15 = high_resolution_clock::now();
        
        for(int i{0};i<length;i++){
            

            
            // 5 random integers to pick the index from the old t values
            //oldTValsIndex = RNG.randInt(0,(length-1),5);
            
            
           //HOW LONG DOES THIS TAKE 
            for(int j{0};j<5;j++){
                
                oldTVals[j] = oldthdist[RNG.randInt(0,length-1)];
                if(USE_THETA==0){
                    oldRVals[j] = sqrt(1-(oldTVals[j] * oldTVals[j]));
                }
                //tvalstest[j] = cos(oldTVals[j]);
                //rvalstest[j] = sqrt(1-(tvalstest[j]*tvalstest[j]));
                //cout << oldTVals[j] << " " << oldRVals[j] << endl;
            }
            
            // generate renormalised t value based on input t values, and other predefined parameters
            
            //thdist[i] = renormalise(angleVector,{oldTVals[(5* i)],oldTVals[(5* i)+1],oldTVals[(5* i)+2],oldTVals[(5* i)+3],oldTVals[(5* i)]+4},RNG.randDouble(0,twopi,8),inputs);
            if(i==-1){
                for(int g{0};g<5;g++){
                    tvalstest[g] = cos(oldTVals[g]);
                    rvalstest[g] = sqrt(1-(tvalstest[g] * tvalstest[g]));
                }
                std::cout << matrixReturnTRI(angleVector,oldTVals,RNG.randDouble(0,twopi,8)) << std::endl;
                auto start = high_resolution_clock::now();
                Matrix<complex<double>,20,20> testsys = matrixReturnTRI(angleVector,oldTVals,RNG.randDouble(0,twopi,8));
                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop - start);
                std::cout << (duration.count()) <<  " matrixReturnTRI"<<  std::endl;
                auto start1 = high_resolution_clock::now();
                Matrix<complex<double>,20,1> testinputs = inputVectorReturnTRI(angleVector, oldTVals, inputs);
                auto stop1 = high_resolution_clock::now();
                auto duration1 = duration_cast<microseconds>(stop1 - start1);
                std::cout << (duration1.count()) <<  " inputVectorReturnTRI"<<  std::endl;
                
                auto start2 = high_resolution_clock::now();
                Matrix<complex<double>,20,1> testinv = testsys.ldlt().solve(testinputs);
                auto stop2 = high_resolution_clock::now();
                auto duration2 = duration_cast<microseconds>(stop2 - start2);
                std::cout << (duration2.count()) <<  " inverting"<<  std::endl;

                auto start3 = high_resolution_clock::now();
                Matrix<complex<double>,10,10> testoriginalsys = matrixReturnORIGINALT(tvalstest,rvalstest,RNG.randDouble(0,twopi,8));
                auto stop3 = high_resolution_clock::now();
                auto duration3 = duration_cast<microseconds>(stop3 - start3);
                std::cout << (duration3.count()) <<  " matrixReturnORIGINALT"<<  std::endl;
                
                auto start4 = high_resolution_clock::now();
                Matrix<complex<double>,10,1> testoriginalinputs = inputVectorReturnORIGINALT(angleVector,oldTVals,inputs);
                auto stop4 = high_resolution_clock::now();
                auto duration4 = duration_cast<microseconds>(stop4 - start4);
                std::cout << (duration4.count()) <<  " inputVectorReturnORIGINALT"<<  std::endl;
                
                auto start5 = high_resolution_clock::now();
                Matrix<complex<double>,10,1> testoriginalinv = testoriginalsys.ldlt().solve(testoriginalinputs);
                auto stop5 = high_resolution_clock::now();
                auto duration5 = duration_cast<microseconds>(stop5 - start5);
                std::cout << (duration5.count()) <<  " inverting"<<  std::endl;
                
                auto start6 = high_resolution_clock::now();
                tdist[i] = renormaliseORIGINALTNOINVERSEFASTERMAYBE(tvalstest,rvalstest,{RNG.randDouble(0,twopi),RNG.randDouble(0,twopi), RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi)},inputs);
                auto stop6 = high_resolution_clock::now();
                auto duration6 = duration_cast<microseconds>(stop6 - start6);
                std::cout << (duration6.count()) <<  " individual random vals ORIGINAL"<<  std::endl;

                auto start7 = high_resolution_clock::now();
                tdist[i] = renormaliseORIGINALTNOINVERSEFASTERMAYBE(tvalstest,rvalstest, RNG.randDouble(0,twopi,8),inputs);
                auto stop7 = high_resolution_clock::now();
                auto duration7 = duration_cast<microseconds>(stop7 - start7);
                std::cout << (duration7.count()) <<  " grouped random vals ORIGINAL"<<  std::endl;

                auto start8 = high_resolution_clock::now();
                thdist[i] = renormaliseTRI(angleVector,oldTVals, RNG.randDouble(0,twopi,8),inputs);
                auto stop8 = high_resolution_clock::now();
                auto duration8 = duration_cast<microseconds>(stop8 - start8);
                std::cout << (duration8.count()) <<  " grouped random vals TRI"<<  std::endl;


                auto start9 = high_resolution_clock::now();
                tdist[i] = renormaliseTRI(angleVector,oldTVals,{RNG.randDouble(0,twopi),RNG.randDouble(0,twopi), RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi),RNG.randDouble(0,twopi)},inputs);
                auto stop9 = high_resolution_clock::now();
                auto duration9 = duration_cast<microseconds>(stop9 - start9);
                std::cout << (duration9.count()) <<  " individual random vals TRI"<<  std::endl;
                



            }
            //HOW LONG DOES THIS TAKE
            //thdist[i] = renormaliseTRI(angleVector,oldTVals,RNG.randDouble(0,twopi,8),inputs);
            //std::cout << thdist[i] << std::endl;
            //HOW LONG DO THESE TAKE
            if(USE_THETA==1){
                //cout << oldTVals[0] << oldTVals[1] << endl;
                thdist[i] = renormaliseTRI(angleVector, oldTVals, RNG.randDouble(0,twopi,8), inputs);
               // cout << thdist[i] << endl;
                tdist[i] = std::cos(thdist[i]);
                gdist[i] = tdist[i] * tdist[i];
                zdist[i] = std::log((1/gdist[i])-1);
                //cout << tdist[i] << " " << zdist[i] << endl;
            } else{
                //cout << renormaliseORIGINALTNOINVERSE(oldTVals, oldRVals, RNG.randDouble(0,twopi,8), inputs) << endl;
                //starts[i] = high_resolution_clock::now();
                tdist[i] = renormaliseORIGINALTNOINVERSE(oldTVals, oldRVals, RNG.randDouble(0,twopi,8), inputs);
                //stops[i] = high_resolution_clock::now();
                //cout << "t renormalisation time: " << duration_cast<microseconds>(stops[i]-starts[i]).count() << endl;
                //cout << oldTVals[0] << "  " << oldRVals[0] << endl;
                //cout << oldTVals[1] << "  " << oldRVals[1] << endl;
                //cout << oldTVals[2] << "  " << oldRVals[2] << endl;
                //cout << oldTVals[3] << "  " << oldRVals[3] << endl;
                //cout << oldTVals[4] << "  " << oldRVals[4] << endl;
                thdist[i] = std::acos(tdist[i]);
                gdist[i] = std::pow(tdist[i],2);
                zdist[i] = std::log((1/gdist[i])-1);
                //cout << tdist[i] << " " << zdist[i]  << endl;
            }
            }

        auto stop15 = high_resolution_clock::now();
        auto duration15 = duration_cast<microseconds>(stop15 - start15);

        std::cout << (duration15.count()) << std::endl;

        if(DEBUG_MODE && current_rank == 0){
            std::cout <<std::endl << "Renormalised! Now creating new distributions from data" << std::endl;
        }
        // create histogram from new data

        int zbinlength = 2 * zbound / zbinsize;
        vector< double> newbinz(zbinlength);
        vector< double> newbinzrev(zbinlength);
        vector< double> newbinth;
        vector< double> newbint;
        vector< double> newbing;
        newbinz = binCounts(zdist,-zbound,zbound,zbinsize, length);
        for(int i{25000};i<25100;i++){
            cout << newbinz[i] << endl;
        }

        if(symmetrise){
            if(DEBUG_MODE && current_rank == 0){
                std::cout << "The z distribution is going to be symmetrised" <<std::endl;
            }
            // Each z distribution value in the first half is taken to be the arithmetic mean of that value and the corresponding value at the other end of the distribution
            std::reverse_copy(std::begin(newbinz),std::end(newbinz),std::begin(newbinzrev));
            vector< double> symdist(length);
            for(int i{0};i< zbinlength;i++){
                newbinz[i] = (newbinz[i] +newbinzrev[i])/2;
            }
            if(DEBUG_MODE && current_rank == 0){
                std::cout << "Laundering symmetrised distribution.." <<std::endl;
            }
            symdist = launder(newbinz,-zbound,zbound,length,zbinsize, RNG);
            if(DEBUG_MODE && current_rank == 0){
                std::cout << "..Done!" <<std::endl;
            }
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

        for(int i{0};i<newbinth.size();i++){
            allbins[i] = newbinth[i];
        }
        for(int i{0};i<newbint.size();i++){
            allbins[i+newbinth.size()] = newbint[i];
        }
        for(int i{0};i<newbing.size();i++){
            allbins[i+newbinth.size()+newbint.size()] = newbing[i];
        }
        for(int i{0};i<newbinz.size();i++){
            allbins[i+newbinth.size()+newbint.size()+newbing.size()] = newbinz[i];
        }


        if(current_rank!= 0){
            std::cout << "Sending bins from rank: " << current_rank << std::endl;
            MPI_Send(&allbins[0],allbins.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
        } else {
            for(int i{0};i<allbins.size();i++){
                avgbins[i] = allbins[i];
            }
            for(int i{1};i<number_of_processes;i++){
                std::cout << "Receiving bins to rank: "<< current_rank << std::endl;
                MPI_Recv(&allbins[0],allbins.size(),MPI_DOUBLE,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
                for(int j{0};j<allbins.size();j++){
                    avgbins[j]+=allbins[j];
                }
            }
            for(int j{0};j<allbins.size();j++){
                allbins[j] = avgbins[j]/number_of_processes;
            }
        }
        if(current_rank==0){

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
                    std::ofstream allbinsout (path + outputPath + "dists" + std::to_string(k+1) + ".txt");
                    allbinsout << binsth.size() << std::endl;
                    allbinsout << binst.size() << std::endl;
                    allbinsout << binsg.size() <<std::endl;
                    allbinsout << binsz.size() << std::endl;
                    for(int i{0};i<binsth.size();i++){
                        binsth[i] = allbins[i];
                        allbinsout << binsth[i] << std::endl;
                    }
                    for(int i{0};i<binst.size();i++){
                        binst[i] = allbins[i+binsth.size()];
                        allbinsout << binst[i] << std::endl;
                    }
                    for(int i{0};i<binsg.size();i++){
                        binsg[i] = allbins[i+binsth.size()+binsg.size()];
                        allbinsout << binsg[i] << std::endl;
                    }
                    for(int i{0};i<binsz.size();i++){
                        binsz[i] = allbins[i+binsth.size()+binsg.size()+binst.size()];
                        allbinsout << binsz[i] << std::endl;
                    }
                    allbinsout.close();
                }
            }
        }
        // std::cout << inv <<std::endl;
        //std::cout << system1 << std::endl;
    }
    MPI_Finalize();
    
    
    
    
    return 0;
}

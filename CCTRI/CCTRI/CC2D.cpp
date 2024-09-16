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
// 0 -> const in th, 1-> linear in t;
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
string getPath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  string fullpath = string( result, (count > 0) ? count : 0 );
  return fullpath.substr(0,fullpath.size()-filename.size());
}


template <typename Out>
void split(const string &s, char delim, Out result) {
    istringstream iss(s);
    string item;
    while (getline(iss, item, delim)) {
        *result++ = item;
    }
}

std::vector<string> split(const string &s, char delim) {
    std::vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}


int main(int argc, char* argv[])
{

    int current_rank;
    int initialisation_error;
    int number_of_processes;
    MPI_Status status;
    //MPI INITIALISATION

    initialisation_error = MPI_Init(&argc,&argv);
    if(initialisation_error != 0){
        cout << "\n";
        cout << "Fatal error!\n";
        cout << "MPI_Init returned ierr = " << initialisation_error << "\n";
        exit(1);
    }

    //Get the number of processes

    MPI_Comm_size(MPI_COMM_WORLD,&number_of_processes);

    //Determine rank of current process

    MPI_Comm_rank(MPI_COMM_WORLD,&current_rank);

    cout << "Initialised process number: " << current_rank << endl;

    //Each process creates their own number generator with varying seed
    mt19937_64 re(seed + current_rank);
    RNG.gen = re;

    /*argument numbers
    noOfSamples -> argv[1] -> arguments[0]
    noOfSteps -> argv[2] -> arguments[1]
    offsetVal -> argv[3] -> arguments[2]
    symmetrise -> argv[4] -> arguments[3]
    readIn -> argv[5] -> arguments[4]
    readInAddress -> argv[6] -> arguments[5]
    "CCTRI-[length-input]-[steps]-[angleInput]/[finalstep]"
    */
    if (argc < 1) {
        // report version
        cout << argv[0] << " Version " << TRIRG_VERSION_MAJOR << "."
              << TRIRG_VERSION_MINOR << endl;
        cout << "Usage: " << argv[0] << " noOfSamples (int)   noOFSteps (int)  offsetVal (double) multiply/divide (bool) spinAngle (as in 2pi/spinAngle, double)    symmetrise (bool)   readIn (bool)   readInAddress (str)" << endl;
        cout << argv[0] << argv[2] << argv[2] << argv[3] << argv[4] << argv[5] << argv[6] << endl;
        return 1;
    }
    // INITIALISATION PARAMS
    // -------------------------------------------
    std::vector<string> arguments = split(argv[1],' ');
    const  double zbound = 25;
    const  double zbinsize = 0.001;
    const  double thgtbinsize = 0.001;
    vector< double> inputs{1,0,0,0};
    // input length is given as an exponent, input 5 will mean the length is 10^5
    const int lengthInput = std::stoi(arguments[0]);
    const int length = pow(10,stoi(arguments[0]));
    int step = 0;
    // number of renormalisation steps is also read in as an argument
    const int steps = stoi(arguments[1]);
    bool symmetrise =stoi(arguments[3]);
    const  double offsetVal = stod(arguments[2]);
    int readIn = stoi(arguments[4]);
    string readInAddress = arguments[5];
    if(DEBUG_MODE && current_rank == 0){
        cout << "Z bound: " << zbound << endl;
        cout << "Z bin size: " << zbinsize << endl;
        cout << "Unit bin size: " << thgtbinsize << endl;
        cout << "Amount of samples: " << length << endl;
        cout << "Amount of steps: "  << steps << endl;
        cout << "Symmetrise set?: " << symmetrise << endl;
        cout << "Offset value: " << offsetVal << endl;
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
        cout << "distribution and histogram arrays initialised, fetching data path" << endl;
    }
    string path = getPath();
    if(DEBUG_MODE && current_rank == 0){
        cout << "file path found!" <<endl;
    }

    // If the user sets readIn to 1, then an ifstream is opened at the readInAddress (which should be a z distribution)
    if(readIn==1){
        if(DEBUG_MODE && current_rank == 0){
            cout<< "read-in acknowledged, opening input stream" <<endl;
        }
        ifstream currzdist;
        cout << "Opening from" << path << "/Data/" << readInAddress << endl;
        currzdist.open(path + "/Data/" + readInAddress);
        string element;
        int z=0;
        int skips = 1;
        cout << "test" << endl;
        for(int i{0};i<3;i++){
          getline(currzdist,element);
          cout << element << endl;
          skips+= stoi(element);

        }
        for(int i{0};i<skips;i++){
          getline(currzdist,element);
        }

        while(currzdist >> element){
            binsz[z++] = stod(element);
        }
        currzdist.close();
        if(DEBUG_MODE && current_rank == 0){
            cout << "histogram file successfully read in, now generating new distribution based on the histogram file" << endl;
        }
        //Launder is invoked to create a set of samples from the distribution data that was read in
        zdist = launder(binsz,-zbound,zbound,length,zbinsize, RNG);
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+exp(zdist[i]));
            tdist[i] = sqrt(gdist[i]);
            thdist[i] = acos(tdist[i]);
        }
        if(DEBUG_MODE && current_rank == 0){
            cout << "Distibutions created successfully" <<endl;
        }
        //Otherwise standard initialised data is produced, set here to be a uniform distribution in theta.
    }else if(readIn==2){
        if(DEBUG_MODE && current_rank == 0){
            cout<< "read-in acknowledged, opening input stream for raw data" <<endl;
        }
        ifstream currzraw;
        currzraw.open(path + "/Data/" + readInAddress + +"/rawfinal" + ".txt");
        string element;
        int i=0;
        while(currzraw >> element){
            zdist[i++] = stod(element);
        }
        currzraw.close();
        if(DEBUG_MODE && current_rank == 0){
            cout << "raw file successfully read in, now converting to other distributions based on the raw file" << endl;
        }
        //Normal conversions between z and th, t and g
        for(int i{0};i<length;i++){
            gdist[i] = 1/(1+exp(zdist[i]));
            tdist[i] = sqrt(gdist[i]);
            thdist[i] = acos(tdist[i]);
        }

    }else{
        if(DEBUG_MODE && current_rank == 0){
            cout << "Creating random distribution from scratch, uniform in theta" <<endl;
        }
        if(INITIAL_DISTRIBUTION==0){
          for(int i{0};i<length;i++){
            thdist[i] = RNG.randDouble(0,twopi/4);
            tdist[i] = cos(thdist[i]);
            gdist[i] = pow(tdist[i],2);
            zdist[i] = log((1/gdist[i])-1);
          }
        } else {
          for(int i{0};i<length;i++){
            tdist[i] = sqrt(RNG.randDouble(0,1));
            thdist[i] = acos(tdist[i]);
            gdist[i] = pow(tdist[i],2);
            zdist[i] = log((1/gdist[i])-1);
        }
        }
        if(DEBUG_MODE && current_rank == 0){
            cout << "Distributions successfully created!" <<endl;
        }
    }

    //**************************************************
    //  LAUNDER TEST
    //**************************************************
    if(DEBUG_MODE && current_rank==0){
        cout << "Testing Launder function" <<endl;
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
        cout << "Test Passed!" <<endl;
    }
    //**************************************************
    //  APPLY OFFSETVAL
    //**************************************************
    if(offsetVal != 0){
        for(int i{0};i<length;i++){
            zdist[i] = zdist[i] + offsetVal;
            gdist[i] = 1/(exp(zdist[i])+1);
            tdist[i] =sqrt(gdist[i]);
            thdist[i] = acos(tdist[i]);
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
        cout << "Counting bins and creating histograms from created distributions.." <<endl;
    }
    binsth = binCounts(thdist,0,twopi/4,thgtbinsize, length);
    binst = binCounts(tdist,0,1,thgtbinsize, length);
    binsg = binCounts(gdist,0,1,thgtbinsize, length);
    binsz = binCounts(zdist,-zbound,zbound,zbinsize,length);
    if(DEBUG_MODE && current_rank==0){
        cout << "..Done!" <<endl <<endl;
    }
    //create directories to save data to
    string outputPath =  "/Data/CCTRI-"+to_string(lengthInput) + "-" + to_string(steps)  + "-" + to_string(1) +   "/" + to_string(offsetVal).substr(to_string(offsetVal).find(".")+1,3) + "/";
    //fs::current_path(fs::temp_directory_path());
    if(DEBUG_MODE && current_rank==0){
        cout << "Creating directories.." <<endl;
    }
    if(current_rank==0){
        fs::create_directories("." + outputPath);
    }
    //for(int i{0};i<steps+1;i++){
    //    fs::create_directory("." + outputPath + "/" + std::to_string(i));
    //}
   if(DEBUG_MODE && current_rank==0){
        cout << "..Done!" <<endl;
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
        cout << "Sending bins from rank: " << current_rank << endl;
        MPI_Send(&allbins[0],allbins.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
    } else {
        for(int i{0};i<allbins.size();i++){
            avgbins[i] = allbins[i];
        }
        for(int i{1};i<number_of_processes;i++){
            cout << "Receiving bins to rank: "<< current_rank << endl;
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
        cout << "Saving to: " + outputPath << endl;
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
            for(int i{0};i<binsz.size();i++){
                binsz[i] = allbins[i+binsth.size()+binst.size()+binsg.size()];
            }
            ofstream outputth (path + outputPath +"thdist0"+ ".txt");
            ofstream outputt (path + outputPath + "tdist0" +  ".txt");
            ofstream outputg (path + outputPath + "gdist0" +  ".txt");
            ofstream outputz (path + outputPath + "zdist0" + ".txt");

            //std::ofstream rawth (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawth" + std::to_string(0) + ".txt");
            //std::ofstream rawt (path +"/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+ "/rawt" + std::to_string(0) + ".txt");
            //std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(0)+"/rawg" + std::to_string(0) + ".txt");
            ofstream rawz (outputPath +"rawz0" +".txt");
            if(DEBUG_MODE && current_rank==0){
                cout << "Writing to files.." <<endl;
            }
            //write to theta file
            for(int i{0};i<binsth.size();i++){
                outputth << binsth[i] << endl;
            //    std::cout << bint[i] << std::endl;
            }
            // write to t and g files
            for(int i{0};i<binst.size();i++){
                outputt << binst[i] << endl;
                outputg << binsg[i] << endl;
            }
            // write to z file
            for(int i{0};i<binsz.size();i++){
                outputz << binsz[i] << endl;
            }
            // close all files
            outputth.close();
            outputt.close();
            outputg.close();
            outputz.close();
            if(WRITE_OUT_RAW==1){
                for(int i{0};i<length;i++){
                rawz << zdist[i] << endl;
                }
            }


            rawz.close();
            if(DEBUG_MODE && current_rank==0){
                cout << "..Done!" <<endl;
            }
        }else{

            ofstream allbinsout (path + outputPath + "dists0.txt");
            allbinsout << binsth.size() << endl;
            allbinsout << binst.size() << endl;
            allbinsout << binsg.size() <<endl;
            allbinsout << binsz.size() << endl;
            for(int i{0};i<allbins.size();i++){
                allbinsout << allbins[i] << endl;
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
        if(DEBUG_MODE && current_rank == 0){
            cout << "Starting " << to_string(k+1) <<"th iteration" <<endl;
        }
        // initialise new array of theta values
        vector< double> oldzdist(length);
        vector< double> oldtdist(length);
        oldzdist = launder(binsz,-zbound,zbound,length,zbinsize, RNG);
        for(int i{0};i<length;i++){
            oldtdist[i] = sqrt(1/(exp(oldzdist[i])+1));
        }
        if(DEBUG_MODE && current_rank == 0){
            cout << "Renormalising" <<endl;
        }
        // create new t values from old values
        //vector<double> oldTVals(%* length);
        //oldTVals = launder(binsth,0,twopi/4,5 * length,thgtbinsize);
        //vector<int> oldTValsIndex(5);
        vector< double> oldTVals(5);
        vector<double> oldRVals(5);

        vector<high_resolution_clock::time_point> starts(length);
        vector<high_resolution_clock::time_point> stops(length);

        auto start = high_resolution_clock::now();
        //cout << matrixReturnTRI(,  vector<double> t, vector<double> x)
        for(int i{0};i<length;i++){




            // 5 random integers to pick the index from the old t values
            //oldTValsIndex = RNG.randInt(0,(length-1),5);
            //HOW LONG DOES THIS TAKE
            for(int j{0};j<5;j++){

                oldTVals[j] = oldtdist[RNG.randInt(0,length-1)];
                oldRVals[j] = sqrt(1-oldTVals[j] * oldTVals[j]);
            }
            // generate renormalised t value based on input t values, and other predefined parameters

            //thdist[i] = renormalise(angleVector,{oldTVals[(5* i)],oldTVals[(5* i)+1],oldTVals[(5* i)+2],oldTVals[(5* i)+3],oldTVals[(5* i)]+4},RNG.randDouble(0,twopi,8),inputs);
            //HOW LONG DOES THIS TAKE
            starts[i] = high_resolution_clock::now();
            tdist[i] = renormaliseORIGINALT(oldTVals,oldRVals,RNG.randDouble(0,twopi,8),inputs);
            stops[i] = high_resolution_clock::now();
            cout << tdist[i] << endl;
            //std::cout << thdist[i] << std::endl;
            //thdist[i] = acos(tdist[i]);
            gdist[i] = tdist[i] * tdist[i];
            zdist[i] = log((1/gdist[i])-1);

            cout << "renormalisation: " <<duration_cast<microseconds>(stops[i]-starts[i]).count() << endl;
            cout << oldTVals[0] << "  " << oldRVals[0] << endl;
            cout << oldTVals[1] << "  " << oldRVals[1] << endl;
            cout << oldTVals[2] << "  " << oldRVals[2] << endl;
            cout << oldTVals[3] << "  " << oldRVals[3] << endl;
            cout << oldTVals[4] << "  " << oldRVals[4] << endl;
            

        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        std::cout << (duration.count()) << std::endl;

        if(DEBUG_MODE && current_rank == 0){
            cout <<endl << "Renormalised! Now creating new distributions from data" << endl;
        }
        // create histogram from new data

        int zbinlength = 2 * zbound / zbinsize;
        vector< double> newbinz(zbinlength);
        vector< double> newbinzrev(zbinlength);
        vector< double> newbinth;
        vector< double> newbint;
        vector< double> newbing;
        newbinz = binCounts(zdist,-zbound,zbound,zbinsize, length);

        if(symmetrise){
            if(DEBUG_MODE && current_rank == 0){
                cout << "The z distribution is going to be symmetrised" <<endl;
            }
            // Each z distribution value in the first half is taken to be the arithmetic mean of that value and the corresponding value at the other end of the distribution
            reverse_copy(begin(newbinz),end(newbinz),begin(newbinzrev));
            vector< double> symdist(length);
            for(int i{0};i< zbinlength;i++){
                newbinz[i] = (newbinz[i] +newbinzrev[i])/2;
            }
            if(DEBUG_MODE && current_rank == 0){
                cout << "Laundering symmetrised distribution.." <<endl;
            }
            symdist = launder(newbinz,-zbound,zbound,length,zbinsize, RNG);
            if(DEBUG_MODE && current_rank == 0){
                cout << "..Done!" <<endl;
            }
            //for(int i{0};i<symdist.size();i++){
            //    std::cout << symdist[i] <<std::endl;
           // }
            for(int i{0};i<length;i++){
                zdist[i] = symdist[i];
                gdist[i] = 1/(exp(zdist[i])+1);
                tdist[i] =sqrt(gdist[i]);
                thdist[i] = acos(tdist[i]);
            }
            newbinth = binCounts(thdist,0,twopi/4,thgtbinsize, length);

            newbint = binCounts(tdist,0,1,thgtbinsize, length);

            newbing = binCounts(gdist,0,1,thgtbinsize, length);
            newbinz = binCounts(zdist,-zbound,zbound,zbinsize,length);


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
            cout << "Sending bins from rank: " << current_rank << endl;
            MPI_Send(&allbins[0],allbins.size(),MPI_DOUBLE,0,1,MPI_COMM_WORLD);
        } else {
            for(int i{0};i<allbins.size();i++){
                avgbins[i] = allbins[i];
            }
            for(int i{1};i<number_of_processes;i++){
                cout << "Receiving bins to rank: "<< current_rank << endl;
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


          for(int i{0};i<binsth.size();i++){
              binsth[i] = allbins[i];
          }
          for(int i{0};i<binst.size();i++){
            binst[i] = allbins[i+binsth.size()];
          }
          for(int i{0};i<binsg.size();i++){
            binsg[i] = allbins[i + binsth.size() + binst.size()];
          }
          for(int i{0};i<binsz.size();i++){
            binsz[i] = allbins[i + binsth.size() + binsg.size() + binsg.size()];
          }
            // open files to write to
            if(((k+1) % OUTPUT_FREQ == 0) || k==0){
                if(SEP_FILE_OUTPUT){
                    cout << "Saving to: " + outputPath << endl;
                    ofstream outputth (path + outputPath +  "thdist" + to_string(k+1) + ".txt");
                    ofstream outputt (path + outputPath +  "tdist" + to_string(k+1) + ".txt");
                    ofstream outputg (path + outputPath +  "gdist" + to_string(k+1) + ".txt");
                    ofstream outputz (path + outputPath +  "zdist" + to_string(k+1) + ".txt");

                    // std::ofstream rawth (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/thraw.txt");
                    // std::ofstream rawt (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/traw.txt");
                    // std::ofstream rawg (path + "/CCTRI-"+std::to_string(lengthInput) + "-" + std::to_string(steps) + "-" + std::to_string((int)angleInput) + "/" + std::to_string(k+1)+ "/graw.txt");
                    ofstream rawz (path + outputPath + "zraw" + to_string(k+1) + ".txt");

                    //write to theta file
                    if(WRITE_OUT_RAW==1){
                        for(int i{0};i<length;i++){
                        //  rawth << thdist[i] << std::endl;
                        //  rawt << tdist[i] << std::endl;
                        //   rawg << gdist[i] << std::endl;
                        rawz << zdist[i] << endl;

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
                        outputth << binsth[i] << endl;
                        //std::cout << std::endl;
                    }
                    // write to t and g files
                    for(int i{0};i<binst.size();i++){
                        outputt << binst[i] << endl;
                        outputg << binsg[i] << endl;
                    }
                    // write to z file
                    for(int i{0};i<binsz.size();i++){
                        outputz << newbinz[i] << endl;
                    }
                    // close all files
                    outputth.close();
                    outputt.close();
                    outputg.close();
                    outputz.close();
                }else{
                    ofstream allbinsout (path + outputPath + "dists" + to_string(k+1) + ".txt");
                    allbinsout << binsth.size() << endl;
                    allbinsout << binst.size() << endl;
                    allbinsout << binsg.size() <<endl;
                    allbinsout << binsz.size() << endl;
                    for(int i{0}; i<allbins.size();i++){
                      allbinsout << allbins[i] << endl;
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

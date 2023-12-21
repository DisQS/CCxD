
#include <random>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>


const double seed = 12;
using std::rand;
using std::mt19937_64;
using std::vector;

mt19937_64 re(seed);



class randNums {
    public:
        // overloaded constructor to allow for taking arrays or just single random numbers
        vector<double> randDouble(double, double, int);
        double randDouble(double,double);
        vector<int> randInt(int,int,int);
        int randInt(int,int);

};

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
vector<int> randInt(int lower, int upper, int length){
    std::uniform_int_distribution<int> unif(lower,upper);
    vector<int> r(length);
    for(int i{0};i<length;i++){
        r[i] = unif(re);
    }
    return r;
}

int randInt(int lower,int upper){
    std::uniform_int_distribution<int> unif(lower,upper);
    int r = unif(re);
    return r;
}
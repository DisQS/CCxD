#include <random>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>


using std::rand;
using std::mt19937_64;
using std::vector;


const double seed = 12;
mt19937_64 re(seed);


/*
---------------------------
randNums
---------------------------
usage:
    class harbouring all random number related functions, with polymorphic methods
*/
class randNums {
    private:
    public:
        // overloaded constructor to allow for taking arrays or just single random numbers
        vector<double> randDouble(double, double, int);
        double randDouble(double,double);
        vector<int> randInt(int,int,int);
        int randInt(int,int);
};


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
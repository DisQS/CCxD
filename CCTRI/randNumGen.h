#include <random>
#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <string>


using std::vector;
using std::mt19937_64;




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
        vector<long double> randDouble(long double, long double, mt19937_64 ,int);
        long double randDouble(long double,long double,mt19937_64 );
        vector<int> randInt(int,int,mt19937_64 ,int);
        int randInt(int,int,mt19937_64 );
};




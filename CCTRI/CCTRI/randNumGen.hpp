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
    std::mt19937_64 gen;
        // overloaded constructor to allow for taking arrays or just single random numbers
        vector< double> randDouble( double,  double, int);
         double randDouble( double, double);
        vector<int> randInt(int,int,int);
        int randInt(int,int);
};




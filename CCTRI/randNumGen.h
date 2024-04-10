#include <random>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>


using std::vector;




/*
---------------------------
randNums
---------------------------
usage:
    class harbouring all random number related functions, with polymorphic methods
*/
#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
class randNums {
    private:
    public:
        // overloaded constructor to allow for taking arrays or just single random numbers
        vector<double> randDouble(double, double, int);
        double randDouble(double,double);
        vector<int> randInt(int,int,int);
        int randInt(int,int);
};
#endif
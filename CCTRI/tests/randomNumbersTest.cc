#include<gtest/gtest.h>

#include"../CCTRI/generalFunctions.cpp"

#include <random>

using namespace std;


using std::mt19937_64;
randNums rngtest;


TEST(RANDOMNUMBERS, INTTESTS){
    mt19937_64 re(5);
    rngtest.gen = re;
    cout << to_string(rngtest.randInt(1,10,1)[0]) << endl;
    EXPECT_EQ(0, rngtest.randInt(1,10,1)[0]);
}


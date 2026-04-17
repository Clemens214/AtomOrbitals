#include <iostream>
#include <cctype>
#include <random>

#include "sample.hpp"

int main()
{   
    int dotCount = 10000;
    int n = 4, l = 2, m = 1;
    Orbital orbit(n, l, m);
    
    const float pMax = 1;
    const float boxHalf = 2;

    std::cout << "Started updating points!" << std::endl;
    std::vector<Point> points = rejectionSample(dotCount, orbit, pMax, boxHalf);
    std::cout << "Finished updating points!" << std::endl;

    return 0;
}
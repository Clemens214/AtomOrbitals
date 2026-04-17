#include <iostream>
#include <cctype>
#include <random>

#include "sample.hpp"

std::vector<Point> computePoints(int n, int l, int m, int dotCount)
{
    Orbital orbit(n, l, m);
    float boxHalf = 5 * n*n;
    float pMax = computePMax(orbit, boxHalf);
    float rMax = computeRMax(orbit, boxHalf, pMax);
    return rejectionSample(dotCount, orbit, pMax, rMax);
}

int main()
{   
    // Initial orbital
    QuantumState qs = { 3, 2, 0};
    Orbital orbit;
    const int dotCount = 5000;

    std::cout << "Started updating points!" << std::endl;
    std::vector<Point> points = computePoints(qs.n, qs.l, qs.m, dotCount);
    std::cout << "Finished updating points!" << std::endl;

    return 0;
}
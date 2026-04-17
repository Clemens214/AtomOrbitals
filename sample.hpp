#ifndef SAMPLE_H
#define SAMPLE_H

#include <random>
#include <iostream>

#include "Orbital.hpp"

struct Point {
    double value = 0;
    double x = 0;
    double y = 0;
    double z = 0;
};

float computePMax(Orbital orbit, float boxHalf, int gridRes = 50)
{
    float pMax = 0;
    float step = (2 * boxHalf) / gridRes;
    for (int ix = 0; ix < gridRes; ++ix) {
        for (int iy = 0; iy < gridRes; ++iy) {
            for (int iz = 0; iz < gridRes; ++iz) {
                float x = -boxHalf + ix * step;
                float y = -boxHalf + iy * step;
                float z = -boxHalf + iz * step;
                pMax = std::max(pMax, (float)orbit.probability(x, y, z));
            }
        }
    }
    return pMax * 1.1;
}

float computeRMax(Orbital orbit, float boxHalf, float pMax, int gridRes = 50)
{
    float RMax = 0;
    float step = (2 * boxHalf) / gridRes;
    for (int ix = 0; ix < gridRes; ++ix) {
        for (int iy = 0; iy < gridRes; ++iy) {
            for (int iz = 0; iz < gridRes; ++iz) {
                float x = -boxHalf + ix * step;
                float y = -boxHalf + iy * step;
                float z = -boxHalf + iz * step;
                if ( 0.05*pMax <= (float)orbit.probability(x, y, z) )
                    RMax = std::max(RMax, (float)sqrt( x*x + y*y + z*z ));
            }
        }
    }
    return RMax * 1.1;
}

// Rejection sampler
std::vector<Point> rejectionSample(int count, Orbital orbit, float pMax, float rMax)
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> spaceDist(-rMax, rMax);
    std::uniform_real_distribution<double> probDist(0.0f, pMax);

    std::vector<Point> points;
    points.reserve(count);
    int pointNum = 0;
    while (points.size() < count) {
        double x = spaceDist(rng);
        double y = spaceDist(rng);
        double z = spaceDist(rng);
        double val = orbit.probability(x, y, z);
        if (probDist(rng) < val) {
            Point point;
            point.value = val;
            point.x = x;
            point.y = y;
            point.z = z;
            points.push_back(point);
            pointNum++;
        }
        // std::cout << "Number of points found: " << pointNum << std::endl;
    }
    return points;
}

#endif // SAMPLE_H
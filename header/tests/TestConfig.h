#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

#include <string>

struct TestConfig {
    std::string datasetPath;
    std::string outputPrefix;
    double beta;
    double gamma;
    double dt;
    int numSteps;
    int minProcs;
    int maxProcs;

    TestConfig(const std::string& path, const std::string& prefix, 
               double b, double g, double deltaT, int steps, int minP, int maxP)
        : datasetPath(path), outputPrefix(prefix), 
          beta(b), gamma(g), dt(deltaT), numSteps(steps), 
          minProcs(minP), maxProcs(maxP) {}
};

#endif // TEST_CONFIG_H
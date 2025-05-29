#ifndef TEST_CONFIG_H
#define TEST_CONFIG_H

#include <string>
#include <vector>

struct TestConfig {
    std::string datasetPath;
    std::string outputPrefix;
    double beta;
    double gamma;
    double dt;
    int numSteps;
    int minProcesses;
    int maxProcesses;

    TestConfig(const std::string& path, const std::string& prefix, 
              double b, double g, double delta_t, int steps, 
              int minP, int maxP)
        : datasetPath(path), outputPrefix(prefix), 
          beta(b), gamma(g), dt(delta_t), numSteps(steps),
          minProcesses(minP), maxProcesses(maxP) {}
};

class TestSuite {
public:
    std::vector<TestConfig> configs;
    
    void addTemporalTests();
    void addParameterSensitivityTests();
    void addAllTests() {
        addTemporalTests();
        addParameterSensitivityTests();
    }
};

#endif
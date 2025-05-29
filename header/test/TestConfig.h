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
    int minProcs;
    int maxProcs;

    TestConfig(const std::string& path, const std::string& prefix, 
               double b, double g, double deltaT, int steps, int minP, int maxP);
};

class TestSuite {
public:
    std::vector<TestConfig> configs;
    void addAllTests();
    void addTemporalTests();
    void addParameterSensitivityTests();
};

#endif // TEST_CONFIG_H
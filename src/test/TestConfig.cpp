#include "../../header/test/TestConfig.h"
#include <filesystem>

TestConfig::TestConfig(const std::string& path, const std::string& prefix, 
                      double b, double g, double delta_t, int steps, 
                      int minP, int maxP)
    : datasetPath(path), outputPrefix(prefix), 
      beta(b), gamma(g), dt(delta_t), numSteps(steps),
      minProcs(minP), maxProcs(maxP) {}

void TestSuite::addAllTests() {
    addTemporalTests();
    addParameterSensitivityTests();
}

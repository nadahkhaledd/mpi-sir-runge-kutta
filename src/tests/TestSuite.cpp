#include "../../header/tests/TestConfig.h"

void TestSuite::addTemporalTests() {
    const std::string dataDir = "./data/test_datasets/";
    const std::string resultDir = "./data/test_results/";
    
    configs.emplace_back(
        dataDir + "sorted_01-01-2021.csv",
        resultDir + "jan_2021",
        0.3, 0.1, 0.2, 100, 2, 8
    );

    configs.emplace_back(
        dataDir + "sorted_02-05-2021.csv",
        resultDir + "feb_2021",
        0.4, 0.15, 0.2, 100, 2, 8
    );
}

void TestSuite::addParameterSensitivityTests() {
    const std::string dataDir = "./data/test_datasets/";
    const std::string resultDir = "./data/test_results/";

    configs.emplace_back(
        dataDir + "sorted_01-01-2021.csv",
        resultDir + "sensitivity_low_beta",
        0.1, 0.1, 0.2, 100, 4, 4
    );

    configs.emplace_back(
        dataDir + "sorted_01-01-2021.csv",
        resultDir + "sensitivity_high_beta",
        0.5, 0.1, 0.2, 100, 4, 4
    );
}
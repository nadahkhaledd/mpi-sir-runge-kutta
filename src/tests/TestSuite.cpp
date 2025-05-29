#include "../../header/tests/TestConfig.h"

void TestSuite::addTemporalTests() {
    const std::string dataDir = "./data/test_datasets/";
    const std::string resultDir = "./data/test_results/";
    
    // First wave test - January 2021
    configs.emplace_back(
        dataDir + "sorted_01-01-2021.csv",
        resultDir + "jan_2021",
        0.3,    // β - transmission rate
        0.1,    // γ - recovery rate
        0.1,    // dt - time step
        100,    // numSteps
        2,      // minProcs
        8       // maxProcs
    );

    // Second wave test - February 2021
    configs.emplace_back(
        dataDir + "sorted_02-05-2021.csv",
        resultDir + "feb_2021",
        0.3,    // Keep same parameters for fair comparison
        0.1,
        0.1,
        100,
        2,
        8
    );
}

void TestSuite::addParameterSensitivityTests() {
    const std::string dataDir = "./data/test_datasets/";
    const std::string resultDir = "./data/test_results/";

    // Use same dataset but vary parameters
    std::string baseDataset = dataDir + "sorted_01-01-2021.csv";

    // Test different transmission rates (β)
    std::vector<double> betaValues = {0.1, 0.3, 0.5};
    for (double beta : betaValues) {
        configs.emplace_back(
            baseDataset,
            resultDir + "sensitivity_beta_" + std::to_string(beta),
            beta,   // Varying β
            0.1,    // Fixed γ
            0.1,    // Fixed dt
            100,    // Fixed steps
            4,      // Fixed processes for comparison
            4
        );
    }
}
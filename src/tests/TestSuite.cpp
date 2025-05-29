#include "../../header/tests/TestConfig.h"

void TestSuite::addTemporalTests() {
    const std::string resultDir = "./data/test_results/";
    
    // Match the main simulation parameters for proper comparison
    const double dt = 0.2;          // Same as main simulation
    const int numSteps = 100;       // Same as main simulation
    const double beta = 0.3;        // Same as main simulation
    const double gamma = 0.1;       // Same as main simulation

    // Get all CSV files from test dataset directory
    std::filesystem::path dataDir("./data/test_datasets/");
    for (const auto& entry : std::filesystem::directory_iterator(dataDir)) {
        if (entry.path().extension() == ".csv") {
            std::string filename = entry.path().filename().string();
            std::string testName = filename.substr(0, filename.find(".csv"));
            
            configs.emplace_back(
                entry.path().string(),
                resultDir + testName,
                beta,    // Use same β as main
                gamma,   // Use same γ as main
                dt,      // Use same dt as main
                numSteps,// Use same number of steps as main
                2,      // minProcs
                8       // maxProcs
            );
        }
    }
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
#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include "../../header/core/SimulationManager.h"
#include "../../header/core/TimingUtils.h"
#include <filesystem>
#include <iostream>
#include <fstream>

void TestRunner::runTest(const TestConfig& config, MPIHandler& mpi) {
    if (mpi.getRank() == 0) {
        std::filesystem::create_directories("./data/test_results");
        std::cout << "\n=== Running Test ===\n"
                  << "Dataset: " << config.datasetPath << "\n"
                  << "Parameters: β=" << config.beta 
                  << ", γ=" << config.gamma << std::endl;
    }

    auto startTime = TimingUtils::startTimer();
    
    // Use SimulationManager to run the simulation
    auto results = SimulationManager::runSimulation(
        mpi,
        config.datasetPath,
        config.beta,
        config.gamma,
        config.dt,
        config.numSteps,
        config.outputPrefix
    );

    // Save test results
    if (mpi.getRank() == 0) {
        std::string resultsPath = config.outputPrefix + 
                                 "_p" + std::to_string(mpi.getSize()) +
                                 "_results.csv";
        std::ofstream outfile(resultsPath);
        outfile << "Time,S,I,R\n";
        int doublesPerStep = 4;
        for (size_t i = 0; i < results.size(); i += doublesPerStep) {
            outfile << results[i] << "," 
                   << results[i + 1] << ","
                   << results[i + 2] << ","
                   << results[i + 3] << "\n";
        }
    }

    TimingUtils::printAndLogTimingSummary("Test_Total", 
        TimingUtils::stopTimer(startTime), mpi.getRank(), mpi.getSize());
}
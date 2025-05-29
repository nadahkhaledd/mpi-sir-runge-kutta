#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include "../../header/core/SIRModel.h"
#include "../../header/core/GridSimulation.h"
#include "../../header/core/CSVParser.h"
#include "../../header/core/TimingUtils.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <filesystem>

void TestRunner::runTest(const TestConfig& config, MPIHandler& mpi) {
    auto startTime = TimingUtils::startTimer();
    
    if (mpi.getRank() == 0) {
        std::filesystem::create_directories("./data/test_results");
        std::cout << "\n=== Running Test ===\n"
                  << "Configuration:\n"
                  << "- Dataset: " << config.datasetPath << "\n"
                  << "Parameters: β=" << config.beta
                  << ", γ=" << config.gamma
                  << ", dt=" << config.dt << "\n"
                  << "- Steps: " << config.numSteps << "\n"
                  << "- Processes: " << mpi.getSize() << "\n";
    }

    auto setupStartTime = TimingUtils::startTimer();
    // Initialize simulation
    SIRModel model(config.beta, config.gamma, config.dt, config.numSteps);
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    // Load data and run simulation
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData(config.datasetPath);
    }

    // Setup and run simulation
    std::map<int, std::list<int>> allBlocks;
    std::unordered_map<int, std::vector<int>> cellNeighborMap, ghostNeighborMap, blockNeighborMap;
    
    simulation.setupSimulation(mpi, fullData, allBlocks, cellNeighborMap, 
                             ghostNeighborMap, blockNeighborMap);
    TimingUtils::printAndLogTimingSummary("Test_Setup", 
        TimingUtils::stopTimer(setupStartTime), mpi.getRank(), mpi.getSize());

    auto simStartTime = TimingUtils::startTimer();
    // Run the simulation
    auto localResults = simulation.runSimulation();
    TimingUtils::printAndLogTimingSummary("Test_Simulation", 
        TimingUtils::stopTimer(simStartTime), mpi.getRank(), mpi.getSize());

    auto gatherStartTime = TimingUtils::startTimer();
    // Gather results
    std::vector<int> recvCounts, displacements;
    auto globalResults = mpi.gatherResults(localResults, recvCounts, displacements);
    TimingUtils::printAndLogTimingSummary("Test_GatherResults", 
        TimingUtils::stopTimer(gatherStartTime), mpi.getRank(), mpi.getSize());

    // Save results with corrected path (remove extra prefix)
    if (mpi.getRank() == 0) {
        std::string resultsPath = config.outputPrefix + 
                                 "_p" + std::to_string(mpi.getSize()) +
                                 "_results.csv";
        std::ofstream outfile(resultsPath);
        outfile << "Time,S,I,R\n";
        int doublesPerStep = 4;
        for (size_t i = 0; i < globalResults.size(); i += doublesPerStep) {
            outfile << globalResults[i] << "," 
                   << globalResults[i + 1] << ","
                   << globalResults[i + 2] << ","
                   << globalResults[i + 3] << "\n";
        }
        
        std::cout << "Test results saved to: " << resultsPath << "\n";
    }

    TimingUtils::printAndLogTimingSummary("Test_Total", 
        TimingUtils::stopTimer(startTime), mpi.getRank(), mpi.getSize());
}
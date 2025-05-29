#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include "../../header/core/SIRModel.h"
#include "../../header/core/GridSimulation.h"
#include "../../header/core/TimingUtils.h"
#include "../../header/core/CSVParser.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

void TestRunner::runTest(const TestConfig& config, MPIHandler& mpi) {
    if (mpi.getRank() == 0) {
        // Create test results directory if it doesn't exist
        std::filesystem::create_directories("./data/test_results");
        
        std::cout << "\n=== Running Test ===\n"
                  << "Configuration:\n"
                  << "- Dataset: " << config.datasetPath << "\n"
                  << "- Beta: " << config.beta << "\n"
                  << "- Gamma: " << config.gamma << "\n"
                  << "- Time step: " << config.dt << "\n"
                  << "- Steps: " << config.numSteps << "\n"
                  << "- Processes: " << mpi.getSize() << "\n";
    }

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

    auto localResults = simulation.runSimulation();

    // Gather results
    std::vector<int> recvCounts, displacements;
    auto globalResults = mpi.gatherResults(localResults, recvCounts, displacements);

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
}
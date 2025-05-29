#include "../../header/core/SimulationManager.h"
#include "../../header/core/CSVParser.h"
#include "../../header/core/TimingUtils.h"
#include <iostream>

std::vector<double> SimulationManager::runSimulation(
    MPIHandler& mpi,
    const std::string& dataPath,
    double beta,
    double gamma,
    double dt,
    int numSteps,
    const std::string& outputPrefix
) {
    // Initialize model and simulation
    SIRModel model(beta, gamma, dt, numSteps);
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    // Load and broadcast data
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData(dataPath);
        if (fullData.empty()) {
            std::cerr << "Failed to load data from: " << dataPath << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        // Print data validation
        std::cout << "Loaded " << fullData.size() << " cells with " 
                  << (fullData.empty() ? 0 : fullData[0].size()) << " values each" << std::endl;
    }

    // Setup and run simulation
    setupAndRunSimulation(mpi, simulation, fullData);

    // Get results only from rank 0's perspective
    std::vector<double> globalResults;
    auto localResults = simulation.runSimulation();

    // Gather results on rank 0 only
    if (mpi.getRank() == 0) {
        // This will hold the final time series
        globalResults.reserve(numSteps * 4); // time, S, I, R for each step
        
        for (const auto& step : localResults) {
            globalResults.push_back(step[0]); // time
            globalResults.push_back(step[1]); // S
            globalResults.push_back(step[2]); // I
            globalResults.push_back(step[3]); // R
        }
    }

    // Only rank 0's results matter for output
    return globalResults;
}

void SimulationManager::setupAndRunSimulation(
    MPIHandler& mpi,
    GridSimulation& simulation,
    std::vector<std::vector<double>>& fullData
) {
    std::map<int, std::list<int>> allBlocks;
    std::unordered_map<int, std::vector<int>> cellNeighborMap, ghostNeighborMap, blockNeighborMap;
    
    simulation.setupSimulation(
        mpi, fullData, allBlocks, cellNeighborMap, 
        ghostNeighborMap, blockNeighborMap
    );
}

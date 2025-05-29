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

    // Load data
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData(dataPath);
        if (fullData.empty()) {
            std::cerr << "Failed to load data from: " << dataPath << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Setup and run
    setupAndRunSimulation(mpi, simulation, fullData);

    // Get results
    auto localResults = simulation.runSimulation();
    std::vector<int> recvCounts, displacements;
    return mpi.gatherResults(localResults, recvCounts, displacements);
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

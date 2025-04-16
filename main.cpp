#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include <iostream>

int main(int argc, char *argv[]) {
    // Initialize MPI
    MPIHandler mpi(argc, argv);
    
    // Create SIR model with parameters
    SIRModel model(0.3, 0.1, 0.1, 1000);
    
    // Load data (only process 0)
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData("initial_conditions.csv");
        std::cout << "Total rows in input dataset: " << fullData.size() << "\n";
    }
    
    // Distribute data among processes
    std::vector<SIRCell> localGrid = mpi.distributeData(fullData);
    
    // Create and run simulation
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());
    simulation.setGrid(localGrid);
    std::vector<std::vector<double>> localResults = simulation.runSimulation();
    
    // Gather and write results
    std::vector<double> globalResults = mpi.gatherResults(localResults);
    mpi.writeResults(globalResults, localResults.size());
    
    return 0;
}
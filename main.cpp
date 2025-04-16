#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include <iostream>
#include <unordered_map>
#include <map>
#include <list>
#include <vector>
#include <string>

int main(int argc, char *argv[]) {

    const int blockSize=4;
    // Initialize MPI
    MPIHandler mpi(argc, argv);

    // Check for input file argument
    if (argc < 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Create SIR model with parameters
    SIRModel model(0.3, 0.1, 0.1, 1000);

    // Load data (only process 0)
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        std::cout << "Total rows in input dataset: " << fullData.size() << "\n";

        // Create cells and blocks using the sorted dataset
        auto cells = GridSimulation::createCellsMap();
        auto blocks = GridSimulation::divideIntoBlocks(cells, blockSize);

        // Debug: Print blocks
        for (const auto& [blockId, cellList] : blocks) {
            std::cout << "Block " << blockId << ": ";
            for (int cell : cellList) {
                std::cout << cell << " ";
            }
            std::cout << "\n";
        }
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
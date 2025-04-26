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

#include <unordered_map>

std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(int rows, int cols) {
    std::unordered_map<int, std::vector<int>> neighbors;
    for (int i = 0; i < rows * cols; ++i) {
        int row = i / cols;
        int col = i % cols;

        std::vector<int> gridNeighbors;

        if (row > 0) gridNeighbors.push_back((row - 1) * cols + col);     // up
        if (row < rows - 1) gridNeighbors.push_back((row + 1) * cols + col); // down
        if (col > 0) gridNeighbors.push_back(i - 1);                       // left
        if (col < cols - 1) gridNeighbors.push_back(i + 1);               // right

        neighbors[i] = gridNeighbors;
    }
    return neighbors;
}


int main(int argc, char *argv[]) {

    const int blockSize=4;
    // Initialize MPI
    MPIHandler mpi(argc, argv);

    // Check for input file argument
    if (argc < 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Create SIR model with parameters
    SIRModel model(0.3, 0.1, 0.2, 100);

    // Load data (only process 0)
    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        std::cout << "Total rows in input dataset: " << fullData.size() << "\n";

        // Create cells and dynamically divide into optimal blocks
        auto cells = GridSimulation::createCellsMap();
        auto blocks = GridSimulation::divideIntoOptimalBlocks(cells, mpi.getSize());

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
    int rows = 8;
    int cols = 8; // So total = 64

    auto neighborMap = build2DGridNeighborMap(rows, cols);
    simulation.setNeighborMap(neighborMap);
    simulation.setGrid(localGrid);
    std::vector<std::vector<double>> localResults = simulation.runSimulation();

    // Gather and write results
    std::vector<double> globalResults = mpi.gatherResults(localResults);
    mpi.writeResults(globalResults, localResults.size());

    return 0;
}
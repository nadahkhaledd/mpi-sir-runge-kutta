#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include <iostream>
#include <unordered_map>
#include <math.h>
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

// Function to calculate grid dimensions and adjust total cells if necessary
std::pair<int, int> calculateGridDimensions(int totalCells, int numBlocks) {
    // Adjust totalCells to ensure it is divisible by the number of blocks
    if (totalCells % numBlocks != 0) {
        int extraCells = numBlocks - (totalCells % numBlocks);
        totalCells += extraCells; // Add dummy cells only if necessary
        std::cout << "Adjusted total cells to " << totalCells << " to ensure even distribution across blocks.\n";
    } else {
        std::cout << "Total cells (" << totalCells << ") is already evenly divisible by the number of blocks (" << numBlocks << ").\n";
    }

    // Calculate grid dimensions to create a rectangular grid with the optimal number of blocks
    // Try to make the grid as square as possible while respecting the number of blocks
    
    // First, try to arrange blocks in a rectangular grid
    int blockRows = static_cast<int>(std::sqrt(numBlocks));
    while (numBlocks % blockRows != 0) {
        blockRows--;
    }
    int blockCols = numBlocks / blockRows;
    
    // Calculate cells per block
    int cellsPerBlock = totalCells / numBlocks;
    
    // Determine how to arrange cells within each block
    // Try to make each block as square as possible
    int cellsPerBlockRow = static_cast<int>(std::sqrt(cellsPerBlock));
    while (cellsPerBlock % cellsPerBlockRow != 0) {
        cellsPerBlockRow--;
    }
    int cellsPerBlockCol = cellsPerBlock / cellsPerBlockRow;
    
    // Calculate final grid dimensions
    int rows = blockRows * cellsPerBlockRow;
    int cols = blockCols * cellsPerBlockCol;
    
    std::cout << "Grid structure: " << blockRows << "x" << blockCols << " blocks, "
              << "each with " << cellsPerBlockRow << "x" << cellsPerBlockCol << " cells\n";
    std::cout << "Final grid dimensions: " << rows << "x" << cols << " (total cells: " << totalCells << ")\n";
    
    return {rows, cols};
}

// Helper function to flatten a 2D vector into a 1D vector
std::vector<double> flatten(const std::vector<std::vector<double>>& data) {
    std::vector<double> flat;
    for (const auto& row : data) {
        flat.insert(flat.end(), row.begin(), row.end());
    }
    return flat;
}

int main(int argc, char *argv[]) {

    //const int blockSize=4;
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

    // Create and run simulation
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    // Distribute data among processes
    std::vector<SIRCell> localGrid = mpi.distributeData(fullData, [&simulation](const std::vector<double>& rowData) {
        return simulation.mapToSIR(rowData);
    });

    // Ensure localGrid is not empty
    if (localGrid.empty()) {
        std::cerr << "Error: Process " << mpi.getRank() << " received no data.\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Dynamically determine the best number of blocks
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

    // Dynamically calculate grid dimensions
    int totalCells = static_cast<int>(localGrid.size()) * mpi.getSize(); // Total cells across all processes
    int numBlocks = static_cast<int>(blocks.size()); // Use the dynamically determined number of blocks
    auto [rows, cols] = calculateGridDimensions(totalCells, numBlocks);

    auto neighborMap = build2DGridNeighborMap(rows, cols);
    simulation.setNeighborMap(neighborMap);
    simulation.setGrid(localGrid);

    // Run the simulation
    std::vector<std::vector<double>> localResults = simulation.runSimulation();

    // Gather and aggregate results
    std::vector<double> globalFlatResults = mpi.gatherResults(localResults);

    if (mpi.getRank() == 0) {
        // Unflatten global results into a 2D vector
        std::vector<std::vector<double>> globalResults;
        int numColumns = 4; // [time, avg_S, avg_I, avg_R]
        for (size_t i = 0; i < globalFlatResults.size(); i += numColumns) {
            globalResults.push_back({
                globalFlatResults[i], 
                globalFlatResults[i + 1], 
                globalFlatResults[i + 2], 
                globalFlatResults[i + 3]
            });
        }

        // Write results
        mpi.writeResults(globalFlatResults, globalResults.size());
    }

    return 0;
}

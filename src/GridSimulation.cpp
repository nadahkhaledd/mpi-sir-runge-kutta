#include "../header/GridSimulation.h"
#include <mpi.h>
#include <unordered_map>
#include <map>
#include <list>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iostream>

GridSimulation::GridSimulation(const SIRModel& m, int mpiRank, int mpiSize) 
    : model(m), rank(mpiRank), size(mpiSize) {}

void GridSimulation::setGrid(const std::vector<SIRCell>& initialGrid) {
    grid = initialGrid;
}

std::vector<SIRCell>& GridSimulation::getGrid() {
    return grid;
}

int GridSimulation::getLocalSize() const {
    return grid.size();
}


void GridSimulation::setNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {

    neighborMap = map;

}

void GridSimulation::updateGrid() {
    std::vector<SIRCell> newGrid = grid;
    for (size_t i = 0; i < grid.size(); ++i) {
        newGrid[i] = model.rk4Step(grid[i]);
    }
    grid = newGrid;
}

void GridSimulation::updateGridNew() {
    std::vector<SIRCell> newGrid = grid;

    for (size_t i = 0; i < grid.size(); ++i) {
        std::vector<SIRCell> neighbors;

        // Get neighbors from map (if exists)
        if (neighborMap.find(i) != neighborMap.end()) {
            const std::vector<int>& neighborIDs = neighborMap[i];
            for (int j : neighborIDs) {
                if (j >= 0 && j < static_cast<int>(grid.size())) {
                    neighbors.push_back(grid[j]);
                }
            }
        }

        // Use model to compute update using neighbors
        newGrid[i] = model.rk4StepWithNeighbors(grid[i], neighbors);
    }

    grid = newGrid;
}

std::map<std::string, int> GridSimulation::createCellsMap() {
    std::map<std::string, int> cells;
    std::ifstream infile("data/sorted_initial_conditions.csv");
    if (!infile) {
        std::cerr << "Error: Could not open sorted_initial_conditions.csv\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string line;
    std::getline(infile, line); // Skip header
    int cellId = 0;

    while (std::getline(infile, line)) {
        std::istringstream ss(line);
        std::string state;
        std::getline(ss, state, ','); // Assuming the first column is the state name
        cells[state] = cellId++;
    }

    return cells;
}

std::map<int, std::list<int>> GridSimulation::divideIntoBlocks(
    const std::map<std::string, int>& cells, int numBlocks) {
    std::map<int, std::list<int>> blocks;

    // Extract and sort cell IDs
    std::vector<int> sortedCellIds;
    for (const auto& [state, cellId] : cells) {
        sortedCellIds.push_back(cellId);
    }
    std::sort(sortedCellIds.begin(), sortedCellIds.end());

    // Calculate block size dynamically
    int blockSize = sortedCellIds.size() / numBlocks;
    int extra = sortedCellIds.size() % numBlocks;

    // Assign sorted cell IDs to blocks
    int startIndex = 0;
    for (int blockId = 0; blockId < numBlocks; ++blockId) {
        int currentBlockSize = blockSize + (blockId < extra ? 1 : 0); // Distribute remainder evenly
        for (int i = 0; i < currentBlockSize; ++i) {
            blocks[blockId].push_back(sortedCellIds[startIndex++]);
        }
    }

    return blocks;
}

std::map<int, std::list<int>> GridSimulation::divideIntoOptimalBlocks(
    const std::map<std::string, int>& cells, int numProcesses) {
    // Dynamically determine the optimal number of blocks
    int numBlocks = std::max(1, numProcesses * 2); // Example heuristic: 2 blocks per process
    return divideIntoBlocks(cells, numBlocks);
}

std::vector<std::vector<double>> GridSimulation::runSimulation() {
    std::vector<std::vector<double>> results; // [time, avg_S, avg_I, avg_R]
    
    for (int step = 0; step < model.getNumSteps(); ++step) {
        // Update grid
        updateGridNew();
        
        // Compute average S, I, R
        double sumS = 0, sumI = 0, sumR = 0;
        for (auto &cell : grid) {
            sumS += cell.getS();
            sumI += cell.getI();
            sumR += cell.getR();
        }
        
        double avgS = sumS / grid.size();
        double avgI = sumI / grid.size();
        double avgR = sumR / grid.size();
        double timeVal = step * model.getDt();
        
        results.push_back({timeVal, avgS, avgI, avgR});
        
        // Synchronize processes
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    return results;
}
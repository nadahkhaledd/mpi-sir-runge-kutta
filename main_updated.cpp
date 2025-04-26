#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include "header/SIRCell.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(
    int rows, int cols,
    const std::unordered_map<int, int>& cellToBlock,
    std::unordered_map<int, std::vector<int>>& ghostNeighbors
) {
    std::unordered_map<int, std::vector<int>> neighbors;

    int totalCells = rows * cols;
    for (int i = 0; i < totalCells; ++i) {
        if (cellToBlock.find(i) == cellToBlock.end()) continue;

        int row = i / cols;
        int col = i % cols;

        std::vector<int> neighborList;
        std::vector<int> directions = {
            (row > 0)        ? (i - cols) : -1,
            (row < rows - 1) ? (i + cols) : -1,
            (col > 0)        ? (i - 1)    : -1,
            (col < cols - 1) ? (i + 1)    : -1
        };

        for (int neighborId : directions) {
            if (neighborId != -1 && cellToBlock.find(neighborId) != cellToBlock.end()) {
                if (cellToBlock.at(i) == cellToBlock.at(neighborId)) {
                    neighborList.push_back(neighborId);
                } else {
                    ghostNeighbors[i].push_back(neighborId);
                }
            }
        }

        neighbors[i] = neighborList;
    }

    return neighbors;
}

int main(int argc, char *argv[]) {

    const int blockSize = 4;

    MPIHandler mpi(argc, argv);

    if (argc < 1) {
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    SIRModel model(0.3, 0.1, 0.2, 100);

    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        std::cout << "Total rows in input dataset: " << fullData.size() << "\n";
    }

    auto cells = GridSimulation::createCellsMap();
    auto blocks = GridSimulation::divideIntoBlocks(cells, blockSize);

    std::unordered_map<int, int> cellToBlock;
    for (const auto& [blockId, cellList] : blocks) {
        for (int cell : cellList) {
            cellToBlock[cell] = blockId;
        }
    }

    std::vector<SIRCell> localGrid = mpi.distributeData(fullData);

    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    int cols = 10;
    int rows = (fullData.size() + cols - 1) / cols;
    std::unordered_map<int, std::vector<int>> ghostNeighbors;
    auto neighborMap = build2DGridNeighborMap(rows, cols, cellToBlock, ghostNeighbors);

    simulation.setNeighborMap(neighborMap);
    simulation.setGhostNeighborMap(ghostNeighbors);
    simulation.setGrid(localGrid);

    std::vector<std::vector<double>> localResults = simulation.runSimulation();

    std::vector<double> globalResults = mpi.gatherResults(localResults);
    mpi.writeResults(globalResults, localResults.size());

    return 0;
}

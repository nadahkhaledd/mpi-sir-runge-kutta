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
    cellIdToLocalIndex.clear();
    for (int i = 0; i < grid.size(); ++i) {
        cellIdToLocalIndex[grid[i].id] = i;
    }
}

void GridSimulation::exchangeGhostCells() {
    std::unordered_map<int, std::vector<double>> sendBuffers;
    std::unordered_map<int, std::vector<double>> recvBuffers;
    std::vector<MPI_Request> requests;

    // 准备要发送的数据 (发送 S, I, R)
    for (const auto& [cellId, neighborList] : ghostNeighborMap) {
        for (int neighborId : neighborList) {
            int targetRank = neighborId % size;
            if (targetRank != rank) {
                sendBuffers[targetRank].push_back(grid[cellIdToLocalIndex[cellId]].getS());
                sendBuffers[targetRank].push_back(grid[cellIdToLocalIndex[cellId]].getI());
                sendBuffers[targetRank].push_back(grid[cellIdToLocalIndex[cellId]].getR());
            }
        }
    }

    // 准备接收buffer
    for (const auto& [cellId, neighborList] : ghostNeighborMap) {
        for (int neighborId : neighborList) {
            int sourceRank = neighborId % size;
            if (sourceRank != rank) {
                recvBuffers[sourceRank].resize(ghostNeighborMap.size() * 3, 0.0); // 每个ghost cell 3个值
            }
        }
    }

    // 发送
    for (auto& [targetRank, buffer] : sendBuffers) {
        MPI_Request req;
        MPI_Isend(buffer.data(), buffer.size(), MPI_DOUBLE, targetRank, 0, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    // 接收
    for (auto& [sourceRank, buffer] : recvBuffers) {
        MPI_Request req;
        MPI_Irecv(buffer.data(), buffer.size(), MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, &req);
        requests.push_back(req);
    }

    // 等待所有通信完成
    MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

    // 更新ghost cell的SIR
    for (auto& [sourceRank, buffer] : recvBuffers) {
        int idx = 0;
        for (const auto& [cellId, neighborList] : ghostNeighborMap) {
            for (int neighborId : neighborList) {
                if (neighborId % size == sourceRank) {
                    if (cellIdToLocalIndex.find(neighborId) != cellIdToLocalIndex.end()) {
                        int localIdx = cellIdToLocalIndex[neighborId];
                        grid[localIdx].setS(buffer[idx]);
                        grid[localIdx].setI(buffer[idx + 1]);
                        grid[localIdx].setR(buffer[idx + 2]);
                        idx += 3;
                    }
                }
            }
        }
    }
}




std::vector<SIRCell>& GridSimulation::getGrid() {
    return grid;
}

int GridSimulation::getLocalSize() const {
    return grid.size();
}

void GridSimulation::setGhostNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {
    ghostNeighborMap = map;
}

void GridSimulation::setNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {

    neighborMap = map;

}

void GridSimulation::updateGrid() {
    exchangeGhostCells();
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

        int currentCellId = grid[i].id;  // 当前 cell 的 id

        if (neighborMap.find(currentCellId) != neighborMap.end()) {
            const std::vector<int>& neighborIDs = neighborMap[currentCellId];
            for (int neighborId : neighborIDs) {
                if (cellIdToLocalIndex.find(neighborId) != cellIdToLocalIndex.end()) {
                    int neighborIdx = cellIdToLocalIndex[neighborId];
                    neighbors.push_back(grid[neighborIdx]);
                }
                // else 是 ghost neighbor，之后可以补 MPI 通信
            }
        }

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
    const std::map<std::string, int>& cells, int blockSize) {
    std::map<int, std::list<int>> blocks;

    // Extract and sort cell IDs
    std::vector<int> sortedCellIds;
    for (const auto& [state, cellId] : cells) {
        sortedCellIds.push_back(cellId);
    }
    std::sort(sortedCellIds.begin(), sortedCellIds.end());

    // Assign sorted cell IDs to blocks
    int blockId = 0;
    for (size_t i = 0; i < sortedCellIds.size(); ++i) {
        if (i > 0 && i % blockSize == 0) {
            blockId++;
        }
        blocks[blockId].push_back(sortedCellIds[i]);
    }

    return blocks;
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
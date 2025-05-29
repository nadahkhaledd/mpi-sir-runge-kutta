#include "../../header/main/GridSimulation.h"
#include "../../header/main/SIRCell.h"
#include "../../header/main/SIRModel.h"
#include "../../header/main/CSVParser.h"
#include "../../header/main/TimingUtils.h"
#include <mpi.h>
#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <array>
#include <cstdio>
#include <memory>
#include <iomanip>

std::string getCurrentDirectory() {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen("pwd", "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    if (!result.empty() && result.back() == '\n') {
        result.pop_back();
    }
    return result;
}

GridSimulation::GridSimulation(const SIRModel& m, int mpiRank, int mpiSize)
    : model(m), rank(mpiRank), size(mpiSize) {}

std::vector<SIRCell>& GridSimulation::getGrid() {
    return grid;
}

int GridSimulation::getLocalSize() const {
    return static_cast<int>(grid.size());
}

void GridSimulation::setGhostNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {
    ghostNeighborMap = map;
}

void GridSimulation::setGrid(const std::vector<SIRCell>& initialGrid) {
    grid = initialGrid;
    std::cerr << "Rank " << rank << " Warning: setGrid called directly. MPI mappings might be invalid." << std::endl;
    globalToLocalCellIndex.clear();
}

std::map<std::string, int> GridSimulation::createCellsMap() {
    std::cout << "Rank " << rank << ": Starting to create cells map from CSV file.\n";
    std::map<std::string, int> cells;
    std::string filePath;
    try {
        std::string currentDir = getCurrentDirectory();
        filePath = currentDir + "/data/sorted_initial_conditions.csv";
    } catch (const std::runtime_error& e) {
        std::cerr << "Rank " << rank << ": Error getting current directory: " << e.what() << ". Cannot create cells map.\n";
        return {};
    }
    std::ifstream infile(filePath);
    if (!infile) {
        std::cerr << "Rank " << rank << ": Error: Could not open file for createCellsMap: " << filePath << "\n";
        return {};
    }
    std::cout << "Rank " << rank << ": Successfully opened file for createCellsMap: " << filePath << "\n";

    std::string line;
    int cellId = 0;
    bool headerFound = false;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (!headerFound && line.find("Province_State") != std::string::npos) {
            headerFound = true;
            continue;
        }
        if (!headerFound) continue;
        std::istringstream ss(line);
        std::string stateOrIdentifier;
        if (std::getline(ss, stateOrIdentifier, ',')) {
            stateOrIdentifier.erase(0, stateOrIdentifier.find_first_not_of(" \t\n\r\f\v"));
            stateOrIdentifier.erase(stateOrIdentifier.find_last_not_of(" \t\n\r\f\v") + 1);
            if (!stateOrIdentifier.empty()) {
                cells[stateOrIdentifier] = cellId++;
            }
        }
    }
    infile.close();
    std::cout << "Rank " << rank << ": Created cells map with " << cells.size() << " entries.\n";
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
    int totalCells = static_cast<int>(cells.size());
    std::cout << "Rank " << rank << ": Finding optimal block distribution for " << totalCells << " cells...\n";

    // Handle edge case: fewer cells than processes
    if (totalCells <= numProcesses) {
        std::map<int, std::list<int>> blocks;
        int blockId = 0;
        for (const auto& [state, cellId] : cells) {
            blocks[blockId++].push_back(cellId);
            if (blockId >= numProcesses) blockId = 0; // Wrap around
        }
        return blocks;
    }

    // Find all divisors of totalCells
    std::vector<int> divisors;
    for (int i = 1; i <= totalCells; ++i) {
        if (totalCells % i == 0) {
            divisors.push_back(i);
        }
    }

    // Score each configuration
    int bestNumBlocks = 1;
    double bestScore = -1.0;

    for (int blocks : divisors) {
        int cellsPerBlock = totalCells / blocks;

        // Calculate how "square" the grid of blocks would be
        int blockRows = static_cast<int>(std::sqrt(blocks));
        while (blocks % blockRows != 0) {
            blockRows--;
        }
        int blockCols = blocks / blockRows;

        // Calculate how "square" each block would be
        int cellRows = static_cast<int>(std::sqrt(cellsPerBlock));
        while (cellsPerBlock % cellRows != 0) {
            cellRows--;
        }
        int cellCols = cellsPerBlock / cellRows;

        // Calculate final grid dimensions
        int totalRows = blockRows * cellRows;
        int totalCols = blockCols * cellCols;

        // Score based on:
        // 1. How close the grid is to being square (ratio of 1.0 is perfect)
        double gridRatio = static_cast<double>(totalRows) / totalCols;
        if (gridRatio > 1.0) gridRatio = 1.0 / gridRatio; // Ensure ratio is <= 1.0

        // 2. How close each block is to being square
        double blockRatio = static_cast<double>(blockRows) / blockCols;
        if (blockRatio > 1.0) blockRatio = 1.0 / blockRatio;

        // 3. How close each cell block is to being square
        double cellRatio = static_cast<double>(cellRows) / cellCols;
        if (cellRatio > 1.0) cellRatio = 1.0 / cellRatio;

        // 4. Balance between blocks and cells per block
        double balanceFactor = static_cast<double>(blocks) / cellsPerBlock;
        if (balanceFactor > 1.0) balanceFactor = 1.0 / balanceFactor;

        // 5. Penalize extreme configurations
        double penalty = 0.0;
        if (blocks < numProcesses || cellsPerBlock < 5) {
            penalty = 0.3; // Apply a penalty for too few blocks or too few cells per block
        }

        // Combine factors into a single score
        double score = gridRatio * 0.25 + blockRatio * 0.25 + cellRatio * 0.2 + balanceFactor * 0.2 - penalty;

        std::cout << "  Option: " << blocks << " blocks with " << cellsPerBlock
                  << " cells each. Grid: " << totalRows << "x" << totalCols
                  << " Score: " << score << "\n";

        // Select the configuration with the highest score
        if (score > bestScore) {
            bestScore = score;
            bestNumBlocks = blocks;
        }
    }

    int cellsPerBlock = totalCells / bestNumBlocks;
    std::cout << "Rank " << rank << ": Optimal distribution: " << bestNumBlocks << " blocks with "
              << cellsPerBlock << " cells each.\n";

    // Divide cells into the best number of blocks
    return divideIntoBlocks(cells, bestNumBlocks);
}

void GridSimulation::setGridFromLocalData(
    const std::map<int, std::list<int>>& localBlocks,
    const std::map<int, std::vector<double>>& localCellData) {

    grid.clear();
    globalToLocalCellIndex.clear();

    std::vector<int> sortedBlockIds;
    for (const auto& pair : localBlocks) sortedBlockIds.push_back(pair.first);
    std::sort(sortedBlockIds.begin(), sortedBlockIds.end());

    std::vector<int> orderedCellIds;
    for (int blockId : sortedBlockIds) {
        const auto& cellList = localBlocks.at(blockId);
        orderedCellIds.insert(orderedCellIds.end(), cellList.begin(), cellList.end());
    }

    grid.reserve(orderedCellIds.size());
    int localIndex = 0;
    for (int cellId : orderedCellIds) {
        auto it = localCellData.find(cellId);
        if (it != localCellData.end()) {
            grid.push_back(CSVParser::mapToSIR(it->second));
            globalToLocalCellIndex[cellId] = localIndex++;
        }
    }
}

void GridSimulation::setBlockInfo(
    const std::map<int, std::list<int>>& localBlocks,
    const std::unordered_map<int, std::vector<int>>& blockNeighbors) {

    localBlockMap = localBlocks;
    blockNeighborMap = blockNeighbors;

    if (!grid.empty()) {
        std::unordered_map<int, int> tempIndexMap;
        std::vector<int> orderedCellIds;
        std::vector<int> sortedBlockIds;
        for (const auto& pair : localBlocks) sortedBlockIds.push_back(pair.first);
        std::sort(sortedBlockIds.begin(), sortedBlockIds.end());
        for (int blockId : sortedBlockIds) {
            const auto& cellList = localBlocks.at(blockId);
            orderedCellIds.insert(orderedCellIds.end(), cellList.begin(), cellList.end());
        }
        for (size_t i = 0; i < orderedCellIds.size(); ++i) {
            tempIndexMap[orderedCellIds[i]] = static_cast<int>(i);
        }
        globalToLocalCellIndex = std::move(tempIndexMap);
    } else {
        globalToLocalCellIndex.clear();
    }
}

void GridSimulation::setCellNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {
    cellNeighborMap = map;
}

void GridSimulation::setBlockToRankMap(const std::unordered_map<int, int>& map) {
    blockToRankMap = map;
}

void GridSimulation::updateGrid() {
    if (grid.empty()) return;
    std::vector<SIRCell> newGrid = grid;
    for (size_t i = 0; i < grid.size(); ++i) {
        newGrid[i] = model.rk4Step(grid[i]);
    }
    grid = std::move(newGrid);
}

void GridSimulation::updateGridNew() {
    if (grid.empty()) return;
    std::vector<SIRCell> newGrid = grid;
    std::vector<int> localIdxToGlobalId(grid.size(), -1);
    for (const auto& pair : globalToLocalCellIndex) {
        if (pair.second >= 0 && static_cast<size_t>(pair.second) < localIdxToGlobalId.size()) {
            localIdxToGlobalId[pair.second] = pair.first;
        }
    }
    for (size_t i = 0; i < grid.size(); ++i) {
        std::vector<SIRCell> neighbors;
        int globalId = localIdxToGlobalId[i];
        if (globalId != -1 && cellNeighborMap.count(globalId)) {
            for (int neighborGlobalID : cellNeighborMap.at(globalId)) {
                if (globalToLocalCellIndex.count(neighborGlobalID)) {
                    neighbors.push_back(grid[globalToLocalCellIndex.at(neighborGlobalID)]);
                }
            }
        }
        newGrid[i] = model.rk4StepWithNeighbors(grid[i], neighbors);
    }

    // Normalize the grid to ensure S + I + R = 1 for each cell
    for (auto& cell : newGrid) {
        double sum = cell.getS() + cell.getI() + cell.getR();
        if (sum > 0) {
            cell.setS(cell.getS() / sum);
            cell.setI(cell.getI() / sum);
            cell.setR(cell.getR() / sum);
        } else {
            cell.setS(1.0);
            cell.setI(0.0);
            cell.setR(0.0);
        }
    }

    grid = std::move(newGrid);
}

// Main MPI Simulation Loop with boundary exchange
std::vector<std::vector<double>> GridSimulation::runSimulation() {
    MPI_Barrier(MPI_COMM_WORLD); // Sync before overall simulation timing
    double totalSimulationWallTime_Start = TimingUtils::startTimer();

    std::vector<std::vector<double>> localResults;
    // Accumulators for total time spent in each phase by THIS RANK over all steps
    double this_rank_total_mpi_prep_time = 0.0;
    double this_rank_total_mpi_comm_time = 0.0;
    double this_rank_total_local_computation_time = 0.0;

    if (grid.empty()) {
        std::cout << "Rank " << rank << ": Grid is empty. Skipping simulation loop, participating in barriers." << std::endl;
        if (size > 1) {
             for (int step_counter = 0; step_counter < model.getNumSteps(); ++step_counter) { // Renamed step
                 MPI_Barrier(MPI_COMM_WORLD);
             }
        }
        MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks sync before this rank might exit early
        double emptyGridDuration = TimingUtils::stopTimer(totalSimulationWallTime_Start);
        // Each rank with an empty grid prints its own (short) duration and logs it if it's rank 0
        TimingUtils::printAndLogIndividualRankTime("runSimulation_TotalWallTime_EmptyGrid", emptyGridDuration, rank);
        return localResults;
    }

    std::vector<int> localIndexToGlobalId(grid.size()); // Default init, will be populated
    for(const auto& entry : globalToLocalCellIndex) { // Using entry as a pair
        if (entry.second >= 0 && static_cast<size_t>(entry.second) < grid.size()) {
            localIndexToGlobalId[entry.second] = entry.first;
        } else {
             std::cerr << "FATAL Error: Rank " << rank << " Invalid local index " << entry.second << " for global ID " << entry.first << " in runSimulation setup." << std::endl; MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    const int DOUBLES_PER_CELL = 3;
    const int MPI_TAG_COUNT = 100;
    const int MPI_TAG_DATA = 200;
    // std::cout << "Rank " << rank << " starting simulation loop for " << model.getNumSteps() << " steps." << std::endl; // Moved to main.cpp

    // Use a different loop variable name to avoid conflict if 'step' is used inside your condensed MPI code
    for (int step_loop_idx = 0; step_loop_idx < model.getNumSteps(); ++step_loop_idx) {
        // --- 1. Boundary Interaction / MPI Prep ---
        double mpiPrepStartTime_step = TimingUtils::startTimer();
        std::map<int, std::set<int>> localIndicesToSendToRank;
        std::map<int, std::set<int>> globalIdsToReceiveFromRank;
        for (size_t localIndex = 0; localIndex < grid.size(); ++localIndex) {
            int globalId = localIndexToGlobalId[localIndex];
            if (globalId == -1 || !cellNeighborMap.count(globalId)) continue;
            const std::vector<int>& neighborGlobalIds = cellNeighborMap.at(globalId);
            for (int neighborGlobalId : neighborGlobalIds) {
                if (globalToLocalCellIndex.find(neighborGlobalId) == globalToLocalCellIndex.end()) { // If remote
                    int owningRank = -1; bool foundOwner = false;
                    // --- YOUR ORIGINAL placeholder for owner lookup ---
                     for(const auto& btor_entry : blockToRankMap) { 
                        int targetRank = btor_entry.second;
                        if (targetRank == rank) continue;
                        
                        auto blockIt = cellToBlockMap.find(neighborGlobalId);
                        if (blockIt != cellToBlockMap.end() && blockIt->second == btor_entry.first) {
                            owningRank = targetRank;
                            foundOwner = true;
                            break;
                        }
                     }
                    // --- END YOUR ORIGINAL placeholder ---
                    if (foundOwner && owningRank != -1) {
                        localIndicesToSendToRank[owningRank].insert(static_cast<int>(localIndex));
                        globalIdsToReceiveFromRank[owningRank].insert(neighborGlobalId);
                    } else if (size > 1 && !foundOwner && owningRank == -1 ) { /* Optional Warning */ }
                }
            }
        }
        this_rank_total_mpi_prep_time += TimingUtils::stopTimer(mpiPrepStartTime_step);


        // --- 2. MPI Communication (Sizes & Data) ---
        double mpiCommStartTime_step = TimingUtils::startTimer();
        std::map<int, int> sendCounts; std::map<int, int> recvCounts;
        std::vector<MPI_Request> countSendRequests, countRecvRequests;
        std::map<int, int> recvCountBuffers;
        for (const auto& pair_gids : globalIdsToReceiveFromRank) { int neighborRank = pair_gids.first; if (neighborRank >= 0 && neighborRank < size && neighborRank != rank) { recvCounts[neighborRank] = 0; recvCountBuffers[neighborRank] = 0; MPI_Request req = MPI_REQUEST_NULL; MPI_Irecv(&recvCountBuffers[neighborRank], 1, MPI_INT, neighborRank, MPI_TAG_COUNT + neighborRank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) countRecvRequests.push_back(req); } }
        std::map<int, int> sendCountBuffers_map;
        for (const auto& pair_lids : localIndicesToSendToRank) { int neighborRank = pair_lids.first; const auto& localIndicesSet = pair_lids.second; if (neighborRank >= 0 && neighborRank < size && neighborRank != rank) { sendCounts[neighborRank] = localIndicesSet.size(); sendCountBuffers_map[neighborRank] = sendCounts[neighborRank]; MPI_Request req = MPI_REQUEST_NULL; MPI_Isend(&sendCountBuffers_map[neighborRank], 1, MPI_INT, neighborRank, MPI_TAG_COUNT + rank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) countSendRequests.push_back(req); } }
        if (!countRecvRequests.empty()) { MPI_Waitall(countRecvRequests.size(), countRecvRequests.data(), MPI_STATUSES_IGNORE); for(auto const& pair_rcvbuf : recvCountBuffers) { recvCounts[pair_rcvbuf.first] = pair_rcvbuf.second; } }
        if (!countSendRequests.empty()) { MPI_Waitall(countSendRequests.size(), countSendRequests.data(), MPI_STATUSES_IGNORE); }
        
        std::vector<MPI_Request> dataSendRequests_inloop, dataRecvRequests_inloop;
        std::map<int, std::vector<double>> sendDataBuffers_inloop;
        std::map<int, std::vector<double>> recvDataBuffers_inloop;
        std::unordered_map<int, SIRCell> ghostCellData_inloop;

        for (const auto& pair_rcounts : recvCounts) { int neighborRank = pair_rcounts.first; int count = pair_rcounts.second; if (count > 0) { int expectedDoubles = count * DOUBLES_PER_CELL; recvDataBuffers_inloop[neighborRank].resize(expectedDoubles); MPI_Request req = MPI_REQUEST_NULL; MPI_Irecv(recvDataBuffers_inloop[neighborRank].data(), expectedDoubles, MPI_DOUBLE, neighborRank, MPI_TAG_DATA + neighborRank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) dataRecvRequests_inloop.push_back(req); } }
        for (const auto& pair_srank : localIndicesToSendToRank) { int neighborRank = pair_srank.first; const auto& localIndicesSet = pair_srank.second; if (sendCounts.count(neighborRank) && sendCounts[neighborRank] > 0) { int countToSend = sendCounts[neighborRank]; int expectedDoubles = countToSend * DOUBLES_PER_CELL; sendDataBuffers_inloop[neighborRank].reserve(expectedDoubles); for (int localIndex_send : localIndicesSet) { if (localIndex_send >= 0 && static_cast<size_t>(localIndex_send) < grid.size()) { const SIRCell& cell_to_send = grid[localIndex_send]; sendDataBuffers_inloop[neighborRank].push_back(cell_to_send.getS()); sendDataBuffers_inloop[neighborRank].push_back(cell_to_send.getI()); sendDataBuffers_inloop[neighborRank].push_back(cell_to_send.getR()); } else { std::cerr << "FATAL Error: Rank " << rank << " Invalid local index " << localIndex_send << " during send prep."<< std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } } if (sendDataBuffers_inloop[neighborRank].size() != static_cast<size_t>(expectedDoubles)) { std::cerr << "FATAL Error: Rank " << rank << " Send buffer size mismatch rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } MPI_Request req = MPI_REQUEST_NULL; MPI_Isend(sendDataBuffers_inloop[neighborRank].data(), expectedDoubles, MPI_DOUBLE, neighborRank, MPI_TAG_DATA + rank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) dataSendRequests_inloop.push_back(req); } }
        if (!dataRecvRequests_inloop.empty()) { MPI_Waitall(dataRecvRequests_inloop.size(), dataRecvRequests_inloop.data(), MPI_STATUSES_IGNORE); for (auto& pair_datarcv : recvDataBuffers_inloop) { int neighborRank = pair_datarcv.first; auto& buffer = pair_datarcv.second; if (!globalIdsToReceiveFromRank.count(neighborRank)) continue; const auto& expectedGlobalIdsSet = globalIdsToReceiveFromRank.at(neighborRank); std::vector<int> expectedGlobalIds(expectedGlobalIdsSet.begin(), expectedGlobalIdsSet.end()); if (buffer.size() != static_cast<size_t>(recvCounts[neighborRank] * DOUBLES_PER_CELL)) { std::cerr << "FATAL Error: Rank " << rank << " Received data size mismatch rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } if (expectedGlobalIds.size() != static_cast<size_t>(recvCounts[neighborRank])) { std::cerr << "FATAL Error: Rank " << rank << " Mismatch #IDs rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } for (size_t i_unpack = 0; i_unpack < expectedGlobalIds.size(); ++i_unpack) { int globalId_unpack = expectedGlobalIds[i_unpack]; size_t bufferIdx = i_unpack * DOUBLES_PER_CELL; if (bufferIdx + 2 < buffer.size()) { double s_val = buffer[bufferIdx + 0]; double i_val = buffer[bufferIdx + 1]; double r_val = buffer[bufferIdx + 2]; ghostCellData_inloop[globalId_unpack] = SIRCell(s_val, i_val, r_val); } else { std::cerr << "FATAL Error: Rank " << rank << " Buffer read index OOB rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } } } }
        if (!dataSendRequests_inloop.empty()) { MPI_Waitall(dataSendRequests_inloop.size(), dataSendRequests_inloop.data(), MPI_STATUSES_IGNORE); }
        this_rank_total_mpi_comm_time += TimingUtils::stopTimer(mpiCommStartTime_step);


        // --- 3. Local Computation ---
        double localCompStartTime_step = TimingUtils::startTimer();
        std::vector<SIRCell> newGrid(grid.size());
        for (size_t localIndex = 0; localIndex < grid.size(); ++localIndex) {
            const SIRCell& currentCell = grid[localIndex];
            int globalId = localIndexToGlobalId[localIndex];
            if (globalId == -1) {
                newGrid[localIndex] = currentCell;
                if (!grid.empty()) std::cerr << "Rank " << rank << " Warning: Invalid globalId (-1) for localIndex " << localIndex << " in computation step " << step_loop_idx << std::endl;
                continue;
            }
            std::vector<SIRCell> neighborsForUpdate;
            if (cellNeighborMap.count(globalId)) {
                const std::vector<int>& neighborGlobalIds = cellNeighborMap.at(globalId);
                neighborsForUpdate.reserve(neighborGlobalIds.size());
                for (int neighborGlobalId : neighborGlobalIds) {
                    auto it_local = globalToLocalCellIndex.find(neighborGlobalId);
                    if (it_local != globalToLocalCellIndex.end()) {
                        neighborsForUpdate.push_back(grid[it_local->second]);
                    } else {
                         auto it_ghost = ghostCellData_inloop.find(neighborGlobalId);
                         if (it_ghost != ghostCellData_inloop.end()) {
                            neighborsForUpdate.push_back(it_ghost->second);
                         }
                    }
                }
            }
            newGrid[localIndex] = model.rk4StepWithNeighbors(currentCell, neighborsForUpdate);
        }
        // Apply normalization as part of the computation step (from your updateGridNew)
        for (auto& cell : newGrid) {
            double sum_sir = cell.getS() + cell.getI() + cell.getR();
            if (sum_sir > 1e-9) {
                cell.setS(cell.getS() / sum_sir);
                cell.setI(cell.getI() / sum_sir);
                cell.setR(cell.getR() / sum_sir);
            } else {
                cell.setS(1.0); cell.setI(0.0); cell.setR(0.0);
            }
        }
        grid = std::move(newGrid);
        this_rank_total_local_computation_time += TimingUtils::stopTimer(localCompStartTime_step);

        // Record local results
        double sumS=0,sumI=0,sumR=0;
        if(!grid.empty()){
            for(const auto&c:grid){sumS+=c.getS();sumI+=c.getI();sumR+=c.getR();}
            sumS/=grid.size();sumI/=grid.size();sumR/=grid.size();
        }
        double timeVal=static_cast<double>(step_loop_idx + 1)*model.getDt();
        localResults.push_back({timeVal,sumS,sumI,sumR});

        if (size > 1) {
             MPI_Barrier(MPI_COMM_WORLD); // End of step synchronization
        }
    } // End simulation loop

    // --- Print AND LOG Accumulated Computation Timings for THIS RANK ---
    TimingUtils::printAndLogIndividualRankTime("Total_MPI_Prep_In_Loop", this_rank_total_mpi_prep_time, rank);
    TimingUtils::printAndLogIndividualRankTime("Total_MPI_Comm_In_Loop", this_rank_total_mpi_comm_time, rank);
    TimingUtils::printAndLogIndividualRankTime("Total_Local_Computation", this_rank_total_local_computation_time, rank);


    // --- Overall Simulation Time (Print summary from Rank 0 AND LOG IT) ---
    double totalSimWallTime_duration_local = TimingUtils::stopTimer(totalSimulationWallTime_Start);
    TimingUtils::printAndLogTimingSummary("runSimulation_TotalWallTime", totalSimWallTime_duration_local, rank, size);

    std::cout << "Rank " << rank << " finished simulation loop." << std::endl;
    return localResults;
}

void GridSimulation::setNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {
    neighborMap = map;
}

std::unordered_map<int, std::vector<int>> GridSimulation::build2DGridNeighborMap(
    int rows, int cols,
    const std::unordered_map<int, int>& cellToBlock,
    std::unordered_map<int, std::vector<int>>& ghostNeighbors) {
    std::unordered_map<int, std::vector<int>> neighbors;
    int totalCells = rows * cols;

    for (int i = 0; i < totalCells; ++i) {
        if (cellToBlock.find(i) == cellToBlock.end()) continue;

        int row = i / cols;
        int col = i % cols;

        std::vector<int> neighborList;
        std::vector<int> directions = {
            (row > 0) ? (i - cols) : -1,
            (row < rows - 1) ? (i + cols) : -1,
            (col > 0) ? (i - 1) : -1,
            (col < cols - 1) ? (i + 1) : -1};

        for (int neighborId : directions) {
            if (neighborId != -1 && cellToBlock.find(neighborId) != cellToBlock.end()) {
                if (cellToBlock.at(i) == cellToBlock.at(neighborId)) {
                    neighborList.push_back(neighborId); // local
                } else {
                    ghostNeighbors[i].push_back(neighborId); // ghost
                    neighborList.push_back(neighborId);      // also keep in main list
                }
            }
        }

        neighbors[i] = neighborList;
    }

    return neighbors;
}

std::unordered_map<int, std::vector<int>> GridSimulation::buildBlockNeighborMap(
    const std::map<int, std::list<int>>& allBlocks,
    const std::unordered_map<int, std::vector<int>>& cellNeighborMap) {
    std::unordered_map<int, std::vector<int>> blockNeighborMap;
    if (allBlocks.empty() || cellNeighborMap.empty()) return blockNeighborMap;

    std::unordered_map<int, int> cellToBlockMap;
    for (const auto& [blockId, cellList] : allBlocks) {
        for (int cellId : cellList) {
            cellToBlockMap[cellId] = blockId;
        }
    }

    for (const auto& [blockId, cellList] : allBlocks) {
        std::unordered_set<int> neighborBlockIdsSet;
        for (int cellId : cellList) {
            auto it_neighbors = cellNeighborMap.find(cellId);
            if (it_neighbors != cellNeighborMap.end()) {
                for (int neighborCellId : it_neighbors->second) {
                    auto it_block = cellToBlockMap.find(neighborCellId);
                    if (it_block != cellToBlockMap.end()) {
                        int neighborBlockId = it_block->second;
                        if (neighborBlockId != blockId) {
                            neighborBlockIdsSet.insert(neighborBlockId);
                        }
                    }
                }
            }
        }
        blockNeighborMap[blockId] = std::vector<int>(neighborBlockIdsSet.begin(), neighborBlockIdsSet.end());
    }

    return blockNeighborMap;
}

void GridSimulation::setupSimulation(
    MPIHandler& mpi,
    const std::vector<std::vector<double>>& fullData,
    std::map<int, std::list<int>>& allBlocks,
    std::unordered_map<int, std::vector<int>>& cellNeighborMap,
    std::unordered_map<int, std::vector<int>>& ghostNeighborMap,
    std::unordered_map<int, std::vector<int>>& blockNeighborMap) {
    if (mpi.getRank() == 0) {
        std::cout << "Rank 0: Setting up simulation...\n";
        auto cells = createCellsMap();
        allBlocks = divideIntoOptimalBlocks(cells, mpi.getSize());

        std::cout << "Rank 0: Divided cells into " << allBlocks.size() << " blocks.\n";

        // Initialize cellToBlockMap in setupSimulation
        if (mpi.getRank() == 0) {
            for (const auto& [blockId, cellList] : allBlocks) {
                for (int cell : cellList) {
                    cellToBlockMap[cell] = blockId;
                }
            }
        }

        std::unordered_map<int, int> cellToBlock;
        for (const auto& [blockId, cellList] : allBlocks) {
            for (int cell : cellList) {
                cellToBlock[cell] = blockId;
            }
        }

        int totalCells = fullData.size();
        int cols = std::ceil(std::sqrt(totalCells));
        int rows = (totalCells + cols - 1) / cols;

        cellNeighborMap = build2DGridNeighborMap(rows, cols, cellToBlock, ghostNeighborMap);
        blockNeighborMap = buildBlockNeighborMap(allBlocks, cellNeighborMap);

        std::cout << "Rank 0: Built neighbor maps.\n";
    }

    allBlocks = mpi.distributeBlocks(allBlocks);
    std::cout << "Rank " << rank << ": Received " << allBlocks.size() << " blocks after distribution.\n";

    auto localCellData = mpi.getDataForLocalBlocks(allBlocks, fullData);
    std::cout << "Rank " << rank << ": Received data for " << localCellData.size() << " cells.\n";

    blockNeighborMap = mpi.broadcastBlockNeighborMap(blockNeighborMap);
    std::cout << "Rank " << rank << ": Received block neighbor map.\n";

    setGridFromLocalData(allBlocks, localCellData);
    std::cout << "Rank " << rank << ": Set grid from local data.\n";

    setBlockInfo(allBlocks, blockNeighborMap);
    setCellNeighborMap(cellNeighborMap);
    setGhostNeighborMap(ghostNeighborMap);
    std::cout << "Rank " << rank << ": Completed simulation setup.\n";
}

std::pair<int, int> GridSimulation::calculateGridDimensions(int totalCells, int numBlocks) {
    if (totalCells % numBlocks != 0) {
        int extraCells = numBlocks - (totalCells % numBlocks);
        totalCells += extraCells;
    }

    int blockRows = static_cast<int>(std::sqrt(numBlocks));
    while (numBlocks % blockRows != 0) {
        blockRows--;
    }
    int blockCols = numBlocks / blockRows;

    int cellsPerBlock = totalCells / numBlocks;
    int cellsPerBlockRow = static_cast<int>(std::sqrt(cellsPerBlock));
    while (cellsPerBlock % cellsPerBlockRow != 0) {
        cellsPerBlockRow--;
    }
    int cellsPerBlockCol = cellsPerBlock / cellsPerBlockRow;

    int rows = blockRows * cellsPerBlockRow;
    int cols = blockCols * cellsPerBlockCol;

    return {rows, cols};
}
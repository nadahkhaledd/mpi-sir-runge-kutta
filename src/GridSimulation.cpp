#include "../header/GridSimulation.h"
#include "../header/SIRCell.h"
#include "../header/SIRModel.h"
#include "../header/CSVParser.h" // Needed for mapToSIR

#include <mpi.h>
#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>    // <<< ADDED for cout, cerr
#include <ostream>     // <<< ADDED for endl
#include <fstream>
#include <sstream>
#include <algorithm>   // <<< ADDED for std::sort
#include <numeric>
#include <stdexcept>   // <<< ADDED for std::out_of_range, runtime_error
#include <array>
#include <cstdio>
#include <memory>

// --- getCurrentDirectory Helper (Keep as is) ---
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

// --- Constructor ---
GridSimulation::GridSimulation(const SIRModel& m, int mpiRank, int mpiSize)
    : model(m), rank(mpiRank), size(mpiSize) {} // Correct Implementation

// --- Getters ---
std::vector<SIRCell>& GridSimulation::getGrid() {
    return grid; // Correct Implementation
}

int GridSimulation::getLocalSize() const {
    return static_cast<int>(grid.size()); // Correct Implementation
}

// --- Basic Setup ---
void GridSimulation::setGrid(const std::vector<SIRCell>& initialGrid) {
    grid = initialGrid;
    std::cerr << "Rank " << rank << " Warning: setGrid called directly. MPI mappings might be invalid." << std::endl;
    globalToLocalCellIndex.clear();
}

// --- Static Helpers ---
std::map<std::string, int> GridSimulation::createCellsMap() {
    // ... Full CORRECTED implementation from previous step should be here ...
     std::map<std::string, int> cells;
     std::string filePath;
     try {
        std::string currentDir = getCurrentDirectory();
        filePath = currentDir + "/data/sorted_initial_conditions.csv";
     } catch (const std::runtime_error& e) {
         std::cerr << "Error getting current directory: " << e.what() << ". Cannot create cells map." << std::endl;
         return {};
     }
     std::ifstream infile(filePath);
     if (!infile) {
         std::cerr << "Error: Could not open file for createCellsMap: " << filePath << "\n";
         return {};
     }
     std::string line;
     int cellId = 0;
     bool headerFound = false;
     while (std::getline(infile, line)) {
         if (line.empty()) continue;
         if (!headerFound && line.find("Province_State") != std::string::npos) {
             headerFound = true;
             std::cout << "createCellsMap: Found header -> " << line << std::endl;
             continue;
         }
         if (!headerFound) {
             std::cout << "createCellsMap: Skipping pre-header line -> " << line << std::endl;
             continue;
         }
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
     std::cout << "createCellsMap: Created map with " << cells.size() << " entries." << std::endl;
     return cells; // Make sure return is here
}

std::map<int, std::list<int>> GridSimulation::divideIntoBlocks(
    const std::map<std::string, int>& cells, int blockSize) {
    // ... Full implementation from previous step should be here ...
    std::map<int, std::list<int>> blocks;
    if (blockSize <= 0) {
        std::cerr << "Error: Block size must be positive (" << blockSize << "). Cannot divide into blocks." << std::endl;
        return blocks;
    }
    if (cells.empty()) {
        std::cerr << "Warning: cells map is empty in divideIntoBlocks." << std::endl;
        return blocks;
    }
    std::vector<int> sortedCellIds;
    sortedCellIds.reserve(cells.size());
    for (const auto& pair : cells) {
        sortedCellIds.push_back(pair.second);
    }
    std::sort(sortedCellIds.begin(), sortedCellIds.end()); // Requires <algorithm>
    int blockId = 0;
    for (size_t i = 0; i < sortedCellIds.size(); ++i) {
        blocks[blockId].push_back(sortedCellIds[i]);
        if ((i + 1) % static_cast<size_t>(blockSize) == 0 && (i + 1) < sortedCellIds.size()) {
            blockId++;
        }
    }
    return blocks; // Make sure return is here
}

// --- Setup Methods for Block Simulation ---

// Sets the local grid from assigned blocks and the LOCAL data subset
void GridSimulation::setGridFromLocalData(
    const std::map<int, std::list<int>>& localBlocks,
    const std::map<int, std::vector<double>>& localCellData) {

    grid.clear();
    globalToLocalCellIndex.clear();

    std::vector<int> sortedBlockIds;
    sortedBlockIds.reserve(localBlocks.size());
    for(const auto& pair : localBlocks) {
        sortedBlockIds.push_back(pair.first);
    }
    // USE std::sort which requires <algorithm>
    std::sort(sortedBlockIds.begin(), sortedBlockIds.end()); // <<< Needs <algorithm>

    std::vector<int> orderedCellIds;
    for (int blockId : sortedBlockIds) {
       try {
            const auto& cellList = localBlocks.at(blockId);
            orderedCellIds.insert(orderedCellIds.end(), cellList.begin(), cellList.end());
       // Catch block requires <stdexcept>
       } catch (const std::out_of_range& oor) { // <<< Needs <stdexcept>
            std::cerr << "Rank " << rank << " Error: Block ID " << blockId << " not found in localBlocks during setGridFromLocalData. what(): " << oor.what() << std::endl; // <<< Needs <iostream>, <ostream>
            MPI_Abort(MPI_COMM_WORLD, 1);
       }
    }

    grid.reserve(orderedCellIds.size());
    int localIndex = 0;
    int missingDataCount = 0;
    for (int cellId : orderedCellIds) {
        auto it = localCellData.find(cellId);
        if (it != localCellData.end()) {
            const std::vector<double>& rowData = it->second;
            try {
                grid.push_back(CSVParser::mapToSIR(rowData));
                globalToLocalCellIndex[cellId] = localIndex++;
            } catch (const std::exception& e) { // <<< Needs <stdexcept> for base class
                 std::cerr << "Rank " << rank << " Error processing cellId " << cellId << " with mapToSIR: " << e.what() << std::endl; // <<< Needs <iostream>, <ostream>
                 missingDataCount++;
            }
        } else {
            std::cerr << "Warning: Rank " << rank << " Data for cellId " << cellId
                      << " (expected in local blocks) not found in received localCellData. Skipping cell." << std::endl; // <<< Needs <iostream>, <ostream>
            missingDataCount++;
        }
    }

    std::cout << "Rank " << rank << " set grid from local data. Grid size: " << grid.size()
              << ". Processed " << localBlocks.size() << " blocks."
              << (missingDataCount > 0 ? (" Skipped " + std::to_string(missingDataCount) + " cells due to missing data.") : "")
              << std::endl; // <<< Needs <iostream>, <ostream>, <string>
}

// Stores block info and verifies/rebuilds the index map
void GridSimulation::setBlockInfo(
    const std::map<int, std::list<int>>& localBlocks,
    const std::unordered_map<int, std::vector<int>>& blockNeighbors) {
    // ... Full implementation from previous step should be here ...
    localBlockMap = localBlocks;
    blockNeighborMap = blockNeighbors;
    if (!grid.empty()) {
        std::unordered_map<int, int> tempIndexMap;
        std::vector<int> orderedCellIds;
        std::vector<int> sortedBlockIds;
        sortedBlockIds.reserve(localBlocks.size());
        for(const auto& pair : localBlocks) { sortedBlockIds.push_back(pair.first); }
        std::sort(sortedBlockIds.begin(), sortedBlockIds.end()); // <<< Needs <algorithm>
        for (int blockId : sortedBlockIds) {
            if(localBlocks.count(blockId)){
                const auto& cellList = localBlocks.at(blockId);
                orderedCellIds.insert(orderedCellIds.end(), cellList.begin(), cellList.end());
            }
        }
        if (orderedCellIds.size() != grid.size()) {
            std::cerr << "Error: Rank " << rank << " Mismatch in setBlockInfo. Cell count from localBlocks ("
                      << orderedCellIds.size() << ") does not match current grid size (" << grid.size()
                      << "). Ensure setGridFromLocalData was called first with the correct blocks/data." << std::endl; // <<< Needs <iostream>, <ostream>
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        for(size_t i = 0; i < orderedCellIds.size(); ++i) {
            tempIndexMap[orderedCellIds[i]] = static_cast<int>(i);
        }
        globalToLocalCellIndex = std::move(tempIndexMap);
    } else if (!localBlocks.empty()) {
        globalToLocalCellIndex.clear();
        std::cout << "Warning: Rank " << rank << " setBlockInfo called, grid is empty. globalToLocalCellIndex cleared." << std::endl; // <<< Needs <iostream>, <ostream>
    } else {
         globalToLocalCellIndex.clear();
    }
    std::cout << "Rank " << rank << " set block info. Local blocks: " << localBlockMap.size()
              << ". Neighbor map size: " << blockNeighborMap.size()
              << ". Index map size: " << globalToLocalCellIndex.size() << std::endl; // <<< Needs <iostream>, <ostream>
}

// Stores the cell neighbor map
void GridSimulation::setCellNeighborMap(const std::unordered_map<int, std::vector<int>>& map) {
    // ... Full implementation from previous step should be here ...
    cellNeighborMap = map;
    std::cout << "Rank " << rank << " set cell neighbor map. Size: " << cellNeighborMap.size() << std::endl; // <<< Needs <iostream>, <ostream>
}

// Stores the block-to-rank map
void GridSimulation::setBlockToRankMap(const std::unordered_map<int, int>& map) {
    // ... Full implementation from previous step should be here ...
    blockToRankMap = map;
     std::cout << "Rank " << rank << " set block-to-rank map. Size: " << blockToRankMap.size() << std::endl; // <<< Needs <iostream>, <ostream>
}

// --- Simulation Logic Implementations ---

// Original updateGrid (no neighbors)
void GridSimulation::updateGrid() {
    if (grid.empty()) return;
    std::vector<SIRCell> newGrid = grid;
    for (size_t i = 0; i < grid.size(); ++i) {
        newGrid[i] = model.rk4Step(grid[i]);
    }
    grid = std::move(newGrid);
}

// Update using LOCAL neighbors only (mainly for testing now)
void GridSimulation::updateGridNew() {
    // ... Full implementation from previous step should be here ...
     if (grid.empty()) return;
     std::vector<SIRCell> newGrid = grid;
     std::vector<int> localIdxToGlobalId(grid.size(), -1);
     for(const auto& pair : globalToLocalCellIndex) {
         if(pair.second >= 0 && static_cast<size_t>(pair.second) < localIdxToGlobalId.size()) {
            localIdxToGlobalId[pair.second] = pair.first;
         }
     }
     for (size_t i = 0; i < grid.size(); ++i) {
         std::vector<SIRCell> neighbors;
         int globalId = localIdxToGlobalId[i];
         if (globalId != -1 && cellNeighborMap.count(globalId)) {
             const std::vector<int>& neighborGlobalIDs = cellNeighborMap.at(globalId);
             neighbors.reserve(neighborGlobalIDs.size());
             for (int neighborGlobalID : neighborGlobalIDs) {
                 if (globalToLocalCellIndex.count(neighborGlobalID)) {
                     neighbors.push_back(grid[globalToLocalCellIndex.at(neighborGlobalID)]);
                 }
             }
         }
         newGrid[i] = model.rk4StepWithNeighbors(grid[i], neighbors);
     }
     grid = std::move(newGrid);
}


// Main MPI Simulation Loop with boundary exchange
std::vector<std::vector<double>> GridSimulation::runSimulation() {
    // ... Full implementation from previous step should be here ...
    // Ensure all std::cout, std::cerr, std::endl uses are valid
    // Ensure the placeholder for findRankOwningCell is handled or replaced
    std::vector<std::vector<double>> localResults;
    if (grid.empty()) {
        std::cout << "Rank " << rank << " has no cells. Entering synchronization loop." << std::endl;
        if (size > 1) {
             for (int step = 0; step < model.getNumSteps(); ++step) {
                 MPI_Barrier(MPI_COMM_WORLD);
             }
        }
        return localResults;
    }
    std::vector<int> localIndexToGlobalId(grid.size());
    for(const auto& [globalId, localIndex] : globalToLocalCellIndex) {
        if (localIndex >= 0 && static_cast<size_t>(localIndex) < grid.size()) {
            localIndexToGlobalId[localIndex] = globalId;
        } else {
             std::cerr << "FATAL Error: Rank " << rank << " Invalid local index " << localIndex << " for global ID " << globalId << std::endl; MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    const int DOUBLES_PER_CELL = 3;
    const int MPI_TAG_COUNT = 100;
    const int MPI_TAG_DATA = 200;
    std::cout << "Rank " << rank << " starting simulation loop for " << model.getNumSteps() << " steps." << std::endl;
    for (int step = 0; step < model.getNumSteps(); ++step) {
        std::map<int, std::set<int>> localIndicesToSendToRank;
        std::map<int, std::set<int>> globalIdsToReceiveFromRank;
        for (size_t localIndex = 0; localIndex < grid.size(); ++localIndex) {
            int globalId = localIndexToGlobalId[localIndex];
            if (!cellNeighborMap.count(globalId)) continue;
            const std::vector<int>& neighborGlobalIds = cellNeighborMap.at(globalId);
            for (int neighborGlobalId : neighborGlobalIds) {
                if (globalToLocalCellIndex.find(neighborGlobalId) == globalToLocalCellIndex.end()) {
                    int owningRank = -1; bool foundOwner = false;
                    // --- REPLACE THIS WITH EFFICIENT LOOKUP ---
                     for(const auto& [b_id, b_rank] : blockToRankMap) { if (b_rank == rank) continue; /* if (isCellInBlock(...)) { owningRank=b_rank; foundOwner=true; break; } */ }
                    // --- END REPLACE ---
                    if (foundOwner && owningRank != -1) {
                        localIndicesToSendToRank[owningRank].insert(static_cast<int>(localIndex));
                        globalIdsToReceiveFromRank[owningRank].insert(neighborGlobalId);
                    } else if (size > 1 && !foundOwner) { /* Optional Warning */ }
                }
            }
        }
        std::map<int, int> sendCounts; std::map<int, int> recvCounts;
        std::vector<MPI_Request> countSendRequests, countRecvRequests;
        std::map<int, int> recvCountBuffers;
        for (const auto& pair : globalIdsToReceiveFromRank) { int neighborRank = pair.first; if (neighborRank >= 0 && neighborRank < size && neighborRank != rank) { recvCounts[neighborRank] = 0; recvCountBuffers[neighborRank] = 0; MPI_Request req = MPI_REQUEST_NULL; MPI_Irecv(&recvCountBuffers[neighborRank], 1, MPI_INT, neighborRank, MPI_TAG_COUNT + neighborRank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) countRecvRequests.push_back(req); } }
        std::map<int, int> sendCountBuffers;
        for (const auto& pair : localIndicesToSendToRank) { int neighborRank = pair.first; const auto& localIndicesSet = pair.second; if (neighborRank >= 0 && neighborRank < size && neighborRank != rank) { sendCounts[neighborRank] = localIndicesSet.size(); sendCountBuffers[neighborRank] = sendCounts[neighborRank]; MPI_Request req = MPI_REQUEST_NULL; MPI_Isend(&sendCountBuffers[neighborRank], 1, MPI_INT, neighborRank, MPI_TAG_COUNT + rank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) countSendRequests.push_back(req); } }
        if (!countRecvRequests.empty()) { MPI_Waitall(countRecvRequests.size(), countRecvRequests.data(), MPI_STATUSES_IGNORE); for(auto const& [nRank, bufferVal] : recvCountBuffers) { recvCounts[nRank] = bufferVal; } }
        if (!countSendRequests.empty()) { MPI_Waitall(countSendRequests.size(), countSendRequests.data(), MPI_STATUSES_IGNORE); }
        std::vector<MPI_Request> dataSendRequests, dataRecvRequests;
        std::map<int, std::vector<double>> sendDataBuffers; std::map<int, std::vector<double>> recvDataBuffers;
        std::unordered_map<int, SIRCell> ghostCellData;
        for (const auto& pair : recvCounts) { int neighborRank = pair.first; int count = pair.second; if (count > 0) { int expectedDoubles = count * DOUBLES_PER_CELL; recvDataBuffers[neighborRank].resize(expectedDoubles); MPI_Request req = MPI_REQUEST_NULL; MPI_Irecv(recvDataBuffers[neighborRank].data(), expectedDoubles, MPI_DOUBLE, neighborRank, MPI_TAG_DATA + neighborRank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) dataRecvRequests.push_back(req); } }
        for (const auto& pair : localIndicesToSendToRank) { int neighborRank = pair.first; const auto& localIndicesSet = pair.second; if (sendCounts.count(neighborRank) && sendCounts[neighborRank] > 0) { int countToSend = sendCounts[neighborRank]; int expectedDoubles = countToSend * DOUBLES_PER_CELL; sendDataBuffers[neighborRank].reserve(expectedDoubles); for (int localIndex : localIndicesSet) { if (localIndex >= 0 && static_cast<size_t>(localIndex) < grid.size()) { const SIRCell& cell = grid[localIndex]; sendDataBuffers[neighborRank].push_back(cell.getS()); sendDataBuffers[neighborRank].push_back(cell.getI()); sendDataBuffers[neighborRank].push_back(cell.getR()); } else { std::cerr << "FATAL Error: Rank " << rank << " Invalid local index " << localIndex << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } } if (sendDataBuffers[neighborRank].size() != static_cast<size_t>(expectedDoubles)) { std::cerr << "FATAL Error: Rank " << rank << " Send buffer size mismatch rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } MPI_Request req = MPI_REQUEST_NULL; MPI_Isend(sendDataBuffers[neighborRank].data(), expectedDoubles, MPI_DOUBLE, neighborRank, MPI_TAG_DATA + rank, MPI_COMM_WORLD, &req); if (req != MPI_REQUEST_NULL) dataSendRequests.push_back(req); } }
        if (!dataRecvRequests.empty()) { MPI_Waitall(dataRecvRequests.size(), dataRecvRequests.data(), MPI_STATUSES_IGNORE); for (auto& pair : recvDataBuffers) { int neighborRank = pair.first; auto& buffer = pair.second; if (!globalIdsToReceiveFromRank.count(neighborRank)) continue; const auto& expectedGlobalIdsSet = globalIdsToReceiveFromRank.at(neighborRank); std::vector<int> expectedGlobalIds(expectedGlobalIdsSet.begin(), expectedGlobalIdsSet.end()); if (buffer.size() != static_cast<size_t>(recvCounts[neighborRank] * DOUBLES_PER_CELL)) { std::cerr << "FATAL Error: Rank " << rank << " Received data size mismatch rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } if (expectedGlobalIds.size() != static_cast<size_t>(recvCounts[neighborRank])) { std::cerr << "FATAL Error: Rank " << rank << " Mismatch #IDs rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } for (size_t i = 0; i < expectedGlobalIds.size(); ++i) { int globalId = expectedGlobalIds[i]; size_t bufferIdx = i * DOUBLES_PER_CELL; if (bufferIdx + 2 < buffer.size()) { double s_val = buffer[bufferIdx + 0]; double i_val = buffer[bufferIdx + 1]; double r_val = buffer[bufferIdx + 2]; ghostCellData[globalId] = SIRCell(s_val, i_val, r_val); } else { std::cerr << "FATAL Error: Rank " << rank << " Buffer read index OOB rank " << neighborRank << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); } } } }
        if (!dataSendRequests.empty()) { MPI_Waitall(dataSendRequests.size(), dataSendRequests.data(), MPI_STATUSES_IGNORE); }
        std::vector<SIRCell> newGrid(grid.size());
        for (size_t localIndex = 0; localIndex < grid.size(); ++localIndex) { const SIRCell& currentCell = grid[localIndex]; int globalId = localIndexToGlobalId[localIndex]; std::vector<SIRCell> neighborsForUpdate; if (cellNeighborMap.count(globalId)) { const std::vector<int>& neighborGlobalIds = cellNeighborMap.at(globalId); neighborsForUpdate.reserve(neighborGlobalIds.size()); for (int neighborGlobalId : neighborGlobalIds) { auto it_local = globalToLocalCellIndex.find(neighborGlobalId); if (it_local != globalToLocalCellIndex.end()) { neighborsForUpdate.push_back(grid[it_local->second]); } else { auto it_ghost = ghostCellData.find(neighborGlobalId); if (it_ghost != ghostCellData.end()) { neighborsForUpdate.push_back(it_ghost->second); } } } } newGrid[localIndex] = model.rk4StepWithNeighbors(currentCell, neighborsForUpdate); }
        grid = std::move(newGrid);
        double sumS = 0, sumI = 0, sumR = 0; if (!grid.empty()) { for (const auto& cell : grid) { sumS += cell.getS(); sumI += cell.getI(); sumR += cell.getR(); } sumS /= grid.size(); sumI /= grid.size(); sumR /= grid.size(); } double timeVal = static_cast<double>(step + 1) * model.getDt(); localResults.push_back({timeVal, sumS, sumI, sumR});
        if (size > 1) { MPI_Barrier(MPI_COMM_WORLD); }
    }
    std::cout << "Rank " << rank << " finished simulation loop." << std::endl;
    return localResults; // Make sure return is here
}
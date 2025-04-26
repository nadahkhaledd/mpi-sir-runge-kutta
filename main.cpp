#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h" // Keep for Rank 0 loading
#include "header/GridSimulation.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <unordered_map> // Keep for neighbor maps
#include <unordered_set> // Needed for buildBlockNeighborMap

// --- Restored Function Implementations ---

// Builds neighbor map for a simple 2D grid ( toroidal or not depending on logic )
std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(int rows, int cols) {
    std::unordered_map<int, std::vector<int>> neighbors;
    int totalCells = rows * cols;
    if (totalCells <= 0) return neighbors; // Handle invalid dimensions

    for (int i = 0; i < totalCells; ++i) {
        int row = i / cols;
        int col = i % cols;

        std::vector<int> currentNeighbors;
        currentNeighbors.reserve(4); // Max 4 neighbors in simple grid

        // Calculate potential neighbor indices (Wrap around for toroidal - adjust if not needed)
        // int up = ((row - 1 + rows) % rows) * cols + col;
        // int down = ((row + 1) % rows) * cols + col;
        // int left = row * cols + ((col - 1 + cols) % cols);
        // int right = row * cols + ((col + 1) % cols);

        // Non-toroidal version (stops at edges)
        if (row > 0) currentNeighbors.push_back((row - 1) * cols + col);     // up
        if (row < rows - 1) currentNeighbors.push_back((row + 1) * cols + col); // down
        if (col > 0) currentNeighbors.push_back(row * cols + (col - 1));     // left
        if (col < cols - 1) currentNeighbors.push_back(row * cols + (col + 1)); // right

        // Add diagonal neighbors if needed (adjust logic)
        // ...

        if (!currentNeighbors.empty()) {
            neighbors[i] = currentNeighbors;
        }
    }
    return neighbors; // Make sure to return the map
}

// Builds a map of adjacent blocks based on cell adjacency
std::unordered_map<int, std::vector<int>> buildBlockNeighborMap(
    const std::map<int, std::list<int>>& allBlocks, // Needs blocks map
    const std::unordered_map<int, std::vector<int>>& cellNeighborMap) { // Needs cell neighbor map

    std::unordered_map<int, std::vector<int>> blockNeighborMap;
    if (allBlocks.empty() || cellNeighborMap.empty()) {
        return blockNeighborMap; // Return empty if input is empty
    }

    // Efficiently map cell ID -> block ID
    std::unordered_map<int, int> cellToBlockMap;
    for (const auto& [blockId, cellList] : allBlocks) {
        for (int cellId : cellList) {
            cellToBlockMap[cellId] = blockId;
        }
    }

    // Determine block neighbors
    for (const auto& [blockId, cellList] : allBlocks) {
        std::unordered_set<int> neighborBlockIdsSet; // Use set to automatically handle duplicates
        for (int cellId : cellList) {
            // Check neighbors of this cell
            auto it_neighbors = cellNeighborMap.find(cellId);
            if (it_neighbors != cellNeighborMap.end()) {
                for (int neighborCellId : it_neighbors->second) {
                    // Find which block the neighbor cell belongs to
                    auto it_block = cellToBlockMap.find(neighborCellId);
                    if (it_block != cellToBlockMap.end()) {
                        int neighborBlockId = it_block->second;
                        // Add if it's a different block
                        if (neighborBlockId != blockId) {
                            neighborBlockIdsSet.insert(neighborBlockId);
                        }
                    }
                }
            }
        }
        // Convert set to vector for the final map
        if (!neighborBlockIdsSet.empty()) {
            blockNeighborMap[blockId] = std::vector<int>(neighborBlockIdsSet.begin(), neighborBlockIdsSet.end());
        }
    }
    return blockNeighborMap; // Make sure to return the map
}

// --- End Restored Functions ---


int main(int argc, char *argv[]) {
    const int blockSize = 4;
    MPIHandler mpi(argc, argv);

    // --- Rank 0 Loads Initial Data and Structures ---
    std::vector<std::vector<double>> fullData;
    std::map<int, std::list<int>> allBlocks;
    std::map<std::string, int> cells;
    std::unordered_map<int, std::vector<int>> blockNeighborMap;

    if (mpi.getRank() == 0) {
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        std::cout << "Rank 0: Total rows loaded from dataset: " << fullData.size() << "\n";
        if (fullData.empty()) {
            std::cerr << "Rank 0 Error: Failed to load initial data. Aborting." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cells = GridSimulation::createCellsMap();
        allBlocks = GridSimulation::divideIntoBlocks(cells, blockSize);

        // Debug print blocks
        for (const auto& [blockId, cellList] : allBlocks) {
            std::cout << "Block " << blockId << ": ";
            for (int cell : cellList) { std::cout << cell << " "; }
            std::cout << "\n";
        }

        int rows = 8; int cols = 8; // Grid dimensions for neighbor calculation
        // Use the restored function with arguments
        auto cellNeighborMapForBlocks = build2DGridNeighborMap(rows, cols);
        // Use the restored function with arguments
        blockNeighborMap = buildBlockNeighborMap(allBlocks, cellNeighborMapForBlocks);

        // Debug print block neighbors
         for (const auto& [blockId, neighborBlocks] : blockNeighborMap) {
             std::cout << "Block " << blockId << " neighbors: ";
             for (int neighborBlock : neighborBlocks) { std::cout << neighborBlock << " "; }
             std::cout << "\n";
         }
    }

    // --- Distribute Block Structure and Necessary Data ---
    std::map<int, std::list<int>> localBlocks = mpi.distributeBlocks(allBlocks);
    std::map<int, std::vector<double>> localCellData = mpi.getDataForLocalBlocks(localBlocks, fullData);
    blockNeighborMap = mpi.broadcastBlockNeighborMap(blockNeighborMap);

    // --- Setup Simulation ---
    SIRModel model(0.3, 0.1, 0.2, 100);
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    int rows = 8; int cols = 8; // Ensure consistency
    auto globalCellNeighborMap = build2DGridNeighborMap(rows, cols);
    simulation.setCellNeighborMap(globalCellNeighborMap);

    // Create and broadcast blockToRankMap
    std::unordered_map<int, int> blockToRankMap;
     if (mpi.getRank() == 0) {
         int totalBlocks = allBlocks.size(); // Use size determined by rank 0
         // --- REMOVED unused variables ---
         // int blocksPerProc = totalBlocks / mpi.getSize();
         // int extraBlocks = totalBlocks % mpi.getSize();
         int blockCounter = 0;
          for (const auto& [blockId, _] : allBlocks) { // Iterate in consistent order if needed
              int targetRank = -1;
              // Recalculate distribution logic to find owner rank for each block
              int tempBlocksPerProc = totalBlocks / mpi.getSize();
              int tempExtraBlocks = totalBlocks % mpi.getSize();
              for(int r=0; r<mpi.getSize(); ++r){
                  int rNumBlocks = (r < tempExtraBlocks) ? (tempBlocksPerProc + 1) : tempBlocksPerProc;
                  int rankStartBlockIndex = (r < tempExtraBlocks) ? (r * (tempBlocksPerProc + 1)) : (r * tempBlocksPerProc + tempExtraBlocks);
                  if(blockCounter >= rankStartBlockIndex && blockCounter < rankStartBlockIndex + rNumBlocks){
                      targetRank = r;
                      break;
                  }
              }
              if (targetRank != -1) {
                  blockToRankMap[blockId] = targetRank;
              } else {
                   std::cerr << "Rank 0 Error: Could not determine owner rank for block " << blockId << " (index " << blockCounter << ")" << std::endl;
              }
              blockCounter++;
          }
     }
     // Broadcast blockToRankMap (using temporary manual broadcast)
      int map_size = (mpi.getRank() == 0) ? blockToRankMap.size() : 0;
      MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
      std::vector<int> map_data(map_size * 2); // key, value pairs
      if (mpi.getRank() == 0) {
          int i = 0;
          for (const auto& pair : blockToRankMap) {
              if (i + 1 < map_size * 2) { // Bounds check
                 map_data[i++] = pair.first;
                 map_data[i++] = pair.second;
              } else {
                 std::cerr << "Rank 0 Error: map_data buffer overflow during Bcast prep." << std::endl;
                 break; // Avoid writing out of bounds
              }
          }
      }
      MPI_Bcast(map_data.data(), map_size * 2, MPI_INT, 0, MPI_COMM_WORLD);
      if (mpi.getRank() != 0) {
          blockToRankMap.clear();
          blockToRankMap.reserve(map_size);
          for (int i = 0; i < map_size * 2; i += 2) {
               if (i + 1 < map_size * 2) { // Bounds check
                  blockToRankMap[map_data[i]] = map_data[i + 1];
               } else {
                  std::cerr << "Rank " << mpi.getRank() << " Error: map_data index out of bounds during Bcast unpack." << std::endl;
                  break; // Avoid reading out of bounds
               }
          }
      }
     simulation.setBlockToRankMap(blockToRankMap);


    simulation.setGridFromLocalData(localBlocks, localCellData);
    simulation.setBlockInfo(localBlocks, blockNeighborMap);

    // --- Run Simulation ---
    std::vector<std::vector<double>> localResults = simulation.runSimulation();

    // --- Gather and Write Results ---
    std::vector<double> globalResults = mpi.gatherResults(localResults);
    int numLocalSteps = localResults.size();
    // --- REMOVED unused variable ---
    // int globalSteps = 0;
    mpi.writeResults(globalResults, numLocalSteps);

    return 0;
}
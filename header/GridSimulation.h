// header/GridSimulation.h

#ifndef GRID_SIMULATION_H
#define GRID_SIMULATION_H

#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <string>
#include "SIRCell.h"  // Ensure this path is correct and SIRCell is defined
#include "SIRModel.h" // Ensure this path is correct and SIRModel is defined
#include "MPIHandler.h" // Include MPIHandler for setupSimulation

class GridSimulation {
public:
    // Constructor
    GridSimulation(const SIRModel& m, int mpiRank, int mpiSize);

    static std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(int rows, int cols);
    static std::pair<int, int> calculateGridDimensions(int totalCells, int numBlocks);

    // --- Getters and Basic Setters ---
    std::vector<SIRCell>& getGrid(); // Get modifiable reference to the local grid
    int getLocalSize() const;        // Get number of cells on this process

    // --- Setup Methods ---
    // <<< ADD THIS DECLARATION BACK >>>
    void setGrid(const std::vector<SIRCell>& initialGrid);
    // <<< END ADDED DECLARATION >>>

    // Sets the grid based on assigned blocks and full dataset
    void setGridFromLocalData(
        const std::map<int, std::list<int>>& localBlocks,
        const std::map<int, std::vector<double>>& localCellData);
    // Stores local block info and the global block neighbor map
    void setBlockInfo(const std::map<int, std::list<int>>& localBlocks,
                      const std::unordered_map<int, std::vector<int>>& blockNeighbors);
    // Stores the global cell-to-cell neighbor map
    void setCellNeighborMap(const std::unordered_map<int, std::vector<int>>& map);
    // Stores the global block-to-rank ownership map
    void setBlockToRankMap(const std::unordered_map<int, int>& map);

    // --- Simulation Logic ---
    void exchangeGhostCells();
    // Original update (no neighbors) - Keep or remove depending on need
    void updateGrid();
    // Update with neighbors (used internally by runSimulation, relies on MPI data)
    void updateGridNew(); // This name might be confusing now, consider renaming or removing if unused externally
    
    void setNeighborMap(const std::unordered_map<int, std::vector<int>>& map);
    void setGhostNeighborMap(const std::unordered_map<int, std::vector<int>>& map);
    
    // Main MPI simulation loop performing communication and updates
    std::vector<std::vector<double>> runSimulation();

    // --- Static Helper Methods ---
    // Creates a map from state names (or identifiers) to cell IDs from a file
    static std::map<std::string, int> createCellsMap();
    static std::map<int, std::list<int>> divideIntoBlocks(
        const std::map<std::string, int>& cells, int blockSize);
    static std::map<int, std::list<int>> divideIntoOptimalBlocks(
        const std::map<std::string, int>& cells, int numProcesses);

    static std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(
        int rows, int cols,
        const std::unordered_map<int, int>& cellToBlock,
        std::unordered_map<int, std::vector<int>>& ghostNeighbors);

    static std::unordered_map<int, std::vector<int>> buildBlockNeighborMap(
        const std::map<int, std::list<int>>& allBlocks,
        const std::unordered_map<int, std::vector<int>>& cellNeighborMap);

    void setupSimulation(
        MPIHandler& mpi,
        const std::vector<std::vector<double>>& fullData,
        std::map<int, std::list<int>>& allBlocks,
        std::unordered_map<int, std::vector<int>>& cellNeighborMap,
        std::unordered_map<int, std::vector<int>>& ghostNeighborMap,
        std::unordered_map<int, std::vector<int>>& blockNeighborMap);

private:
    // --- Core Simulation Components ---
    const SIRModel& model;         // Reference to the SIR model parameters and functions
    std::vector<SIRCell> grid;     // Local grid data (SIR state) for cells owned by this process
    int rank;                      // MPI rank of this process
    int size;                      // Total number of MPI processes

    // --- Data Structures for Block-Based MPI Simulation ---
    // Maps local block IDs to the list of global cell IDs within that block
    std::map<int, std::list<int>> localBlockMap;
    // Maps global block IDs to their neighboring global block IDs (full map)
    std::unordered_map<int, std::vector<int>> blockNeighborMap;
    // Maps global cell IDs owned by this process to their index in the local 'grid' vector
    std::unordered_map<int, int> globalToLocalCellIndex;

    // --- Member Variables (Single Declarations) ---
    // Maps global cell IDs to their neighboring global cell IDs (full map)
    std::unordered_map<int, std::vector<int>> cellNeighborMap;
    // Maps global block IDs to the MPI rank that owns them (full map)
    std::unordered_map<int, int> blockToRankMap;
    std::unordered_map<int, std::vector<int>> ghostNeighborMap;
    std::unordered_map<int, int> cellIdToLocalIndex; // global cell ID -> local index
    std::unordered_map<int, std::vector<int>> neighborMap; // local neighbor map if different from cellNeighborMap



}; // End of class GridSimulation definition

#endif // GRID_SIMULATION_H
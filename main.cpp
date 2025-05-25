#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include "header/TimingUtils.h" // <<< INCLUDE TIMING UTILS
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set> // In case your buildXXXMap functions use it
#include <cmath>         // For std::ceil, std::sqrt if buildXXXMap uses it

// Your build2DGridNeighborMap and buildBlockNeighborMap functions should be defined here
// or included if they are in a separate utility file.
// For this example, I'm assuming they are defined in main.cpp or another linked file.
// If not, you'll need to ensure their definitions are available.
// std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(int rows, int cols);
// std::unordered_map<int, std::vector<int>> buildBlockNeighborMap(
// const std::map<int, std::list<int>>& allBlocks,
// const std::unordered_map<int, std::vector<int>>& cellNeighborMap);


int main(int argc, char *argv[])
{
    MPIHandler mpi(argc, argv); // MPI_Init happens here
    // std::cout << "Rank " << mpi.getRank() << ": MPI initialized.\n"; // Moved below TimingUtils init for cleaner log start

    // --- Initialize Timing Log File (Rank 0 only) ---
    if (mpi.getRank() == 0) {
        if (!TimingUtils::initLogFile("./data/timing_log.csv")) { // Or your preferred filename
            std::cerr << "Rank 0: Failed to initialize timing log file. Timing will not be written to file." << std::endl;
            // You might choose to MPI_Abort or continue without file logging
        }
    }
    MPI_Barrier(MPI_COMM_WORLD); // Ensure log file is open before other ranks proceed with timed ops

    std::cout << "Rank " << mpi.getRank() << ": MPI initialized (after timing init).\n";


    std::vector<std::vector<double>> fullData; // Only Rank 0 will populate this meaningfully
    if (mpi.getRank() == 0)
    {
        std::cout << "Rank 0: Loading initial condition data...\n";
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv"); // Ensure path is correct
        if (fullData.empty())
        {
            std::cerr << "Rank 0: Error: Failed to load initial condition data. Aborting.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        std::cout << "Rank 0: Loaded " << fullData.size() << " rows of initial condition data.\n";
    }

    // These maps will be populated by Rank 0 within setupSimulation and then
    // relevant parts distributed or broadcasted.
    std::map<int, std::list<int>> allBlocks_global_on_rank0; // Rank 0 will fill this
    std::unordered_map<int, std::vector<int>> cellCellNeigborMap_global; // Rank 0 will fill this
    std::unordered_map<int, std::vector<int>> ghostCellMap_global;      // Rank 0 will fill this (or it's built by all)
    std::unordered_map<int, std::vector<int>> blockNeighborMap_global;  // Rank 0 will fill this

    // --- Block to Rank Map ---
    // This map needs to be created by Rank 0 based on how blocks are distributed,
    // and then broadcasted to all. The logic for this is typically tied to
    // how distributeBlocks assigns blocks.
    std::unordered_map<int, int> blockToRankMap;
    // (The logic for populating and broadcasting blockToRankMap needs to be here or in setupSimulation)
    // For now, assuming it's handled correctly by setupSimulation or manually broadcasted after.


    SIRModel model(0.3, 0.1, 0.2, 100); // Example parameters
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    // setupSimulation will internally handle:
    // - Rank 0 creating allBlocks_global_on_rank0, cellCellNeigborMap_global, etc.
    // - Calling mpi.distributeBlocks (which will modify allBlocks_global_on_rank0 to be local for each rank)
    // - Calling mpi.getDataForLocalBlocks
    // - Calling mpi.broadcastBlockNeighborMap
    // - Setting up the simulation object's internal state (grid, localBlockMap, etc.)
    std::cout << "Rank " << mpi.getRank() << ": Calling simulation.setupSimulation().\n";
    simulation.setupSimulation(
        mpi,
        fullData,                   // Passed to setupSimulation, used by rank 0
        allBlocks_global_on_rank0,  // Rank 0 populates, mpi.distributeBlocks makes it local
        cellCellNeigborMap_global,  // Rank 0 populates, needs to be broadcasted or built by all
        ghostCellMap_global,        // Rank 0 populates, needs to be broadcasted or built by all
        blockNeighborMap_global     // Rank 0 populates, mpi.broadcastBlockNeighborMap distributes
    );
    std::cout << "Rank " << mpi.getRank() << ": Returned from simulation.setupSimulation().\n";

    // After setupSimulation, allBlocks_global_on_rank0 now holds the *local* blocks for each rank.
    // The simulation object also has its own internal copy of local blocks.

    // !!! REMOVE THE REDUNDANT CALL to mpi.distributeBlocks !!!
    // std::map<int, std::list<int>> localBlocks = mpi.distributeBlocks(allBlocks_global_on_rank0); // This was the issue

    // The blockToRankMap needs to be set. If setupSimulation doesn't handle its broadcast,
    // you need to do it here. Assuming setupSimulation ensures sim object gets it.
    // If not:
    // if (mpi.getRank() == 0) { /* build blockToRankMap based on allBlocks_global_on_rank0 *before* it was made local */ }
    // mpi.broadcastAndSet_blockToRankMap(blockToRankMap, simulation); // Hypothetical function
    simulation.setBlockToRankMap(blockToRankMap); // Make sure blockToRankMap is correctly populated on all ranks


    std::cout << "Rank " << mpi.getRank() << ": Starting simulation...\n";
    std::vector<std::vector<double>> localResults = simulation.runSimulation();
    std::cout << "Rank " << mpi.getRank() << ": Simulation completed.\n";

    std::vector<int> recvCounts, displacements; // Will be populated by gatherResults on Rank 0
    std::vector<double> globalResults = mpi.gatherResults(localResults, recvCounts, displacements);

    mpi.writeResults(globalResults, recvCounts, displacements); // Pass necessary info

    // --- Close Timing Log File (Rank 0 only) ---
    if (mpi.getRank() == 0) {
        TimingUtils::closeLogFile();
    }
    MPI_Barrier(MPI_COMM_WORLD); // Ensure log file is closed before anyone finalizes

    std::cout << "Rank " << mpi.getRank() << ": Finalizing MPI.\n";
    // MPIHandler destructor will call MPI_Finalize()
    return 0;
}
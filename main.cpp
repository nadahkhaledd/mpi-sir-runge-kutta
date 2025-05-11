#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

int main(int argc, char *argv[])
{
    MPIHandler mpi(argc, argv);
    std::cout << "Rank " << mpi.getRank() << ": MPI initialized.\n";

    std::vector<std::vector<double>> fullData;
    if (mpi.getRank() == 0)
    {
        std::cout << "Rank 0: Loading initial condition data...\n";
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        if (fullData.empty())
        {
            std::cerr << "Rank 0: Error: Failed to load initial condition data.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        std::cout << "Rank 0: Loaded " << fullData.size() << " rows of initial condition data.\n";
    }

    std::map<int, std::list<int>> allBlocks;
    std::unordered_map<int, std::vector<int>> cellNeighborMap, ghostNeighborMap, blockNeighborMap;

    SIRModel model(0.3, 0.1, 0.2, 100);
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    simulation.setupSimulation(mpi, fullData, allBlocks, cellNeighborMap, ghostNeighborMap, blockNeighborMap);

    std::map<int, std::list<int>> localBlocks = mpi.distributeBlocks(allBlocks);

    std::cout << "Rank " << mpi.getRank() << ": Starting simulation...\n";
    std::vector<std::vector<double>> localResults = simulation.runSimulation();
    std::cout << "Rank " << mpi.getRank() << ": Simulation completed.\n";

    std::vector<int> recvCounts, displacements;
    std::vector<double> globalResults = mpi.gatherResults(localResults, recvCounts, displacements);

    mpi.writeResults(globalResults, recvCounts, displacements);

    std::cout << "Rank " << mpi.getRank() << ": Finalizing MPI.\n";
    return 0;
}

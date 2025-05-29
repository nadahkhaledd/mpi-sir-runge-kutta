#include <mpi.h>
#include "header/core/MPIHandler.h"
#include "header/core/SIRModel.h"
#include "header/core/CSVParser.h"
#include "header/core/GridSimulation.h"
#include "header/core/TimingUtils.h"
#include "header/core/SimulationManager.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set> 
#include <cmath>
#include <fstream> // For file output

int main(int argc, char *argv[]) {
    MPIHandler mpi(argc, argv);

    // --- Initialize Timing Log File (Rank 0 only) ---
    if (mpi.getRank() == 0) {
        if (!TimingUtils::initLogFile("./data/output/timing_log.csv")) { // Or your preferred filename
            std::cerr << "Rank 0: Failed to initialize timing log file. Timing will not be written to file." << std::endl;
            // You might choose to MPI_Abort or continue without file logging
        }
    }
    MPI_Barrier(MPI_COMM_WORLD); // Ensure log file is open before other ranks proceed with timed ops

    std::cout << "Rank " << mpi.getRank() << ": MPI initialized (after timing init).\n";


    // Run simulation using manager
    auto results = SimulationManager::runSimulation(
        mpi,
        "./data/sorted_initial_conditions.csv",
        0.3, 0.1, 0.2, 100,
        "./data/simulation"
    );

    // Write results
    if (mpi.getRank() == 0) {
        std::ofstream outfile("./data/output/simulation_results.csv");
        outfile << "Time,S,I,R\n";
        int doublesPerStep = 4;
        for (size_t i = 0; i < results.size(); i += doublesPerStep) {
            outfile << results[i] << "," 
                   << results[i + 1] << ","
                   << results[i + 2] << ","
                   << results[i + 3] << "\n";
        }
        TimingUtils::closeLogFile();
    }

    std::cout << "Rank " << mpi.getRank() << ": Finalizing MPI.\n";
    // MPIHandler destructor will call MPI_Finalize()
    return 0;
}
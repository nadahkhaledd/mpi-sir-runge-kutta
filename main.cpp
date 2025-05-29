#include <mpi.h>
#include "header/core/MPIHandler.h"
#include "header/core/SimulationManager.h"
#include "header/core/TimingUtils.h"
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
    MPIHandler mpi(argc, argv);

    if (mpi.getRank() == 0) {
        if (!TimingUtils::initLogFile("./data/output/timing_log.csv")) {
            std::cerr << "Rank 0: Failed to initialize timing log file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto results = SimulationManager::runSimulation(
        mpi,
        "./data/sorted_initial_conditions.csv",
        0.3, 0.1, 0.2, 100,
        "./data/simulation"
    );

    if (mpi.getRank() == 0) {
        std::ofstream outfile("./data/output/simulation_results.csv");
        outfile << "Time,S_avg,I_avg,R_avg\n";  // Simplified header
        int doublesPerStep = 4;
        for (size_t i = 0; i < results.size(); i += doublesPerStep) {
            outfile << results[i] << ","     // Time 
                   << results[i + 1] << ","  // S
                   << results[i + 2] << ","  // I
                   << results[i + 3] << "\n"; // R
        }
        TimingUtils::closeLogFile();
    }

    return 0;
}
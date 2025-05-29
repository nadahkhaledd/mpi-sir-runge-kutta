#include "../../header/test//TestConfig.h"
#include "../../header/test//TestRunner.h"
#include "../../header/main/MPIHandler.h"
#include "../../header/main/SimulationManager.h"
#include "../../header/main/TimingUtils.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <iomanip> // For std::setprecision

void TestRunner::runTest(const TestConfig& config, MPIHandler& mpi) {
    if (mpi.getRank() == 0) {
        std::filesystem::create_directories("./data/test_results");
        std::cout << "\n=== Running Test ===\n"
                  << "Dataset: " << config.datasetPath << "\n"
                  << "Parameters: β=" << config.beta 
                  << ", γ=" << config.gamma 
                  << ", dt=" << config.dt
                  << ", steps=" << config.numSteps << std::endl;
    }

    auto startTime = TimingUtils::startTimer();
    
    // Run simulation and get results
    auto results = SimulationManager::runSimulation(
        mpi,
        config.datasetPath,
        config.beta,
        config.gamma,
        config.dt,
        config.numSteps,
        config.outputPrefix
    );

    // Save test results (only from rank 0)
    if (mpi.getRank() == 0) {
        std::string resultsPath = config.outputPrefix + 
                                "_p" + std::to_string(mpi.getSize()) +
                                "_results.csv";
        
        std::ofstream outfile(resultsPath);
        if (!outfile) {
            std::cerr << "Error: Could not open file " << resultsPath << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        outfile << "Time,S_avg,I_avg,R_avg\n";
        const int valuesPerStep = 4;
        
        for (size_t i = 0; i < results.size(); i += valuesPerStep) {
            if (i + 3 >= results.size()) break;
            outfile << std::fixed << std::setprecision(6)
                   << results[i] << ","
                   << results[i + 1] << ","
                   << results[i + 2] << ","
                   << results[i + 3] << "\n";
        }
        outfile.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    TimingUtils::printAndLogTimingSummary("Test_Total", 
        TimingUtils::stopTimer(startTime), mpi.getRank(), mpi.getSize());
}
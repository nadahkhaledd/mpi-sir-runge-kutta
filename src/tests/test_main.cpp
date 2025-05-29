#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include <mpi.h>  // Add MPI header
#include <iostream>

int main(int argc, char *argv[]) {
    MPIHandler mpi(argc, argv);

    // Initialize test suite with all configurations
    TestSuite testSuite;
    testSuite.addAllTests();
    
    if (mpi.getRank() == 0) {
        std::cout << "=== Starting SIR Model Test Suite ===\n";
        std::cout << "Number of test configurations: " << testSuite.configs.size() << "\n";
    }
    
    // Run each test configuration
    for (const auto& config : testSuite.configs) {
        if (mpi.getRank() == 0) {
            std::cout << "\n=== Test Configuration ===\n" 
                     << "Dataset: " << config.datasetPath << "\n"
                     << "Output: " << config.outputPrefix << "\n"
                     << "Parameters: β=" << config.beta 
                     << ", γ=" << config.gamma
                     << ", dt=" << config.dt
                     << ", steps=" << config.numSteps << "\n"
                     << "Processes: " << config.minProcesses 
                     << "-" << config.maxProcesses << "\n";
        }
        
        TestRunner::runTest(config, mpi);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi.getRank() == 0) {
        std::cout << "\n=== Test Suite Completed ===\n";
    }

    return 0;
}
#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include "../../header/core/TimingUtils.h"
#include <mpi.h>
#include <iostream>
#include <filesystem>

int main(int argc, char *argv[]) {
    MPIHandler mpi(argc, argv);

    if (mpi.getRank() == 0) {
        std::filesystem::create_directories("./data/test_results");
        if (!TimingUtils::initLogFile("./data/test_results/timing_log.csv")) {
            std::cerr << "Failed to initialize timing log file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        std::cout << "\n=== SIR Model Test Suite ===\n";
    }
    
    TestSuite testSuite;
    testSuite.addAllTests();
    
    if (mpi.getRank() == 0) {
        std::cout << "Total test configurations: " << testSuite.configs.size() << "\n\n";
    }
    
    // Run each test configuration
    for (size_t i = 0; i < testSuite.configs.size(); ++i) {
        const auto& config = testSuite.configs[i];
        if (mpi.getRank() == 0) {
            std::cout << "Running test " << (i + 1) << "/" << testSuite.configs.size() 
                     << ": " << config.outputPrefix << std::endl;
        }
        
        TestRunner::runTest(config, mpi);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi.getRank() == 0) {
        TimingUtils::closeLogFile();
        std::cout << "\n=== Test Suite Completed ===\n";
    }

    return 0;
}
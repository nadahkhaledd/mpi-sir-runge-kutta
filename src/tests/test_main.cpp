#include "../../header/tests/TestConfig.h"
#include "../../header/tests/TestRunner.h"
#include "../../header/core/MPIHandler.h"
#include "../../header/core/TimingUtils.h"
#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[]) {
    MPIHandler mpi(argc, argv);

    if (mpi.getRank() == 0) {
        if (!TimingUtils::initLogFile("./data/test_results/timing_log.csv")) {
            std::cerr << "Failed to initialize timing log file." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    TestSuite testSuite;
    testSuite.addAllTests();
    
    if (mpi.getRank() == 0) {
        std::cout << "=== SIR Model Test Suite ===\n";
        std::cout << "Configurations: " << testSuite.configs.size() << "\n";
    }
    
    for (const auto& config : testSuite.configs) {
        TestRunner::runTest(config, mpi);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi.getRank() == 0) {
        TimingUtils::closeLogFile();
        std::cout << "\n=== Test Suite Completed ===\n";
    }

    return 0;
}
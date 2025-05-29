#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

#include "../main/MPIHandler.h"
#include "TestConfig.h"

class TestRunner {
public:
    // Declaration only - implementation moved to cpp file
    static void runTest(const TestConfig& config, MPIHandler& mpi);
};

#endif

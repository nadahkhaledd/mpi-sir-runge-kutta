#ifndef TEST_RUNNER_H
#define TEST_RUNNER_H

#include "../main/MPIHandler.h"
#include "TestConfig.h"

class TestRunner {
public:
    static void runTest(const TestConfig& config, MPIHandler& mpi);
};

#endif

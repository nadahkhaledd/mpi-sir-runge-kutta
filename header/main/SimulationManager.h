#ifndef SIMULATION_MANAGER_H
#define SIMULATION_MANAGER_H

#include "MPIHandler.h"
#include "SIRModel.h"
#include "GridSimulation.h"
#include <string>
#include <vector>

class SimulationManager {
public:
    static std::vector<double> runSimulation(
        MPIHandler& mpi,
        const std::string& dataPath,
        double beta,
        double gamma,
        double dt,
        int numSteps,
        const std::string& outputPrefix
    );

private:
    static void setupAndRunSimulation(
        MPIHandler& mpi,
        GridSimulation& simulation,
        std::vector<std::vector<double>>& fullData
    );
};

#endif

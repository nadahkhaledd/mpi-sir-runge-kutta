#ifndef SIRMODEL_H
#define SIRMODEL_H

#include "SIRCell.h"
#include <vector>

class SIRModel {
private:
    double beta;      // Infection rate
    double gammaRate; // Recovery rate
    double dt;        // Time step
    int numSteps;     // Total simulation steps

public:
    SIRModel(double b = 0.3, double g = 0.1, double timeStep = 0.1, int steps = 1000);
    
    // Getters
    double getBeta() const;
    double getGamma() const;
    double getDt() const;
    int getNumSteps() const;
    
    // RK4 step function for SIR cell dynamics
    SIRCell rk4Step(const SIRCell &current) const;
    SIRCell rk4StepWithNeighbors(const SIRCell& current, const std::vector<SIRCell>& neighbors) const;

};

#endif // SIRMODEL_H
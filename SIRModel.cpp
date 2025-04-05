#include "SIRModel.h"

SIRModel::SIRModel(double b, double g, double timeStep, int steps)
    : beta(b), gammaRate(g), dt(timeStep), numSteps(steps) {}

double SIRModel::getBeta() const { 
    return beta; 
}

double SIRModel::getGamma() const { 
    return gammaRate; 
}

double SIRModel::getDt() const { 
    return dt; 
}

int SIRModel::getNumSteps() const { 
    return numSteps; 
}

SIRCell SIRModel::rk4Step(const SIRCell &current) const {
    double S = current.getS();
    double I = current.getI();
    double R = current.getR();
    
    auto fS = [&](double s, double i) -> double {
        return -beta * s * i;
    };
    auto fI = [&](double s, double i) -> double {
        return beta * s * i - gammaRate * i;
    };
    auto fR = [&](double i) -> double {
        return gammaRate * i;
    };
    
    // K1
    double k1_S = dt * fS(S, I);
    double k1_I = dt * fI(S, I);
    double k1_R = dt * fR(I);
    
    // K2
    double k2_S = dt * fS(S + 0.5 * k1_S, I + 0.5 * k1_I);
    double k2_I = dt * fI(S + 0.5 * k1_S, I + 0.5 * k1_I);
    double k2_R = dt * fR(I + 0.5 * k1_I);
    
    // K3
    double k3_S = dt * fS(S + 0.5 * k2_S, I + 0.5 * k2_I);
    double k3_I = dt * fI(S + 0.5 * k2_S, I + 0.5 * k2_I);
    double k3_R = dt * fR(I + 0.5 * k2_I);
    
    // K4
    double k4_S = dt * fS(S + k3_S, I + k3_I);
    double k4_I = dt * fI(S + k3_S, I + k3_I);
    double k4_R = dt * fR(I + k3_I);
    
    double newS = S + (k1_S + 2*k2_S + 2*k3_S + k4_S) / 6.0;
    double newI = I + (k1_I + 2*k2_I + 2*k3_I + k4_I) / 6.0;
    double newR = R + (k1_R + 2*k2_R + 2*k3_R + k4_R) / 6.0;
    
    return SIRCell(newS, newI, newR);
}
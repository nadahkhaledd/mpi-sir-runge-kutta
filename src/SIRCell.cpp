#include "../header/SIRCell.h"
#include <algorithm>
#include <cmath>

SIRCell::SIRCell(double s, double i, double r) 
    : S(s), I(i), R(r) {
    normalize();
}

double SIRCell::getS() const { 
    return S; 
}

double SIRCell::getI() const { 
    return I; 
}

double SIRCell::getR() const { 
    return R; 
}

void SIRCell::setS(double s) { 
    S = std::max(0.0, std::min(1.0, s)); 
}

void SIRCell::setI(double i) { 
    I = std::max(0.0, std::min(1.0, i)); 
}

void SIRCell::setR(double r) { 
    R = std::max(0.0, std::min(1.0, r)); 
}

void SIRCell::normalize() {
    double sum = S + I + R;
    if (sum > 0) {
        S /= sum;
        I /= sum;
        R /= sum;
    } else {
        // Fallback if all zeros
        S = 0.99;
        I = 0.01;
        R = 0.0;
    }
    
    // Ensure values are valid
    S = std::max(0.0, std::min(1.0, S));
    I = std::max(0.0, std::min(1.0, I));
    R = std::max(0.0, std::min(1.0, R));
}
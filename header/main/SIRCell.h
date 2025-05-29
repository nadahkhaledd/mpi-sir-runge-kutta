#ifndef SIRCELL_H
#define SIRCELL_H

class SIRCell {
private:
    double S, I, R;

public:
    SIRCell(double s = 0.0, double i = 0.0, double r = 0.0);

    // Getters
    double getS() const;
    double getI() const;
    double getR() const;

    // Setters with validation
    void setS(double s);
    void setI(double i);
    void setR(double r);

    // Normalize to ensure S + I + R = 1
    void normalize();
};

#endif // SIRCELL_H
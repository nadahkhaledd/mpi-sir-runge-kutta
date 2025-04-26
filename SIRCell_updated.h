#ifndef SIRCELL_H
#define SIRCELL_H

class SIRCell {
private:
    double S, I, R;

public:
    int id;         
    int x, y;
    int blockID;

    SIRCell(double s = 0.0, double i = 0.0, double r = 0.0, int id_ = -1)
        : S(s), I(i), R(r), id(id_), x(0), y(0), blockID(-1) {}

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
#ifndef GRIDSIMULATION_H
#define GRIDSIMULATION_H

#include <vector>
#include "SIRCell.h"
#include "SIRModel.h"

class GridSimulation {
private:
    std::vector<SIRCell> grid;
    SIRModel model;
    int rank, size;
    
public:
    GridSimulation(const SIRModel& m, int mpiRank, int mpiSize);
    
    void setGrid(const std::vector<SIRCell>& initialGrid);
    
    std::vector<SIRCell>& getGrid();
    
    int getLocalSize() const;
    
    void updateGrid();
    
    std::vector<std::vector<double>> runSimulation();
};

#endif // GRIDSIMULATION_H
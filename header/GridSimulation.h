#ifndef GRIDSIMULATION_H
#define GRIDSIMULATION_H

#include <vector>
#include <unordered_map>
#include <map>
#include <list>
#include <string>
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

    static std::map<std::string, int> createCellsMap();
    static std::map<int, std::list<int>> divideIntoBlocks(
        const std::map<std::string, int>& cells, int blockSize);
};

#endif // GRIDSIMULATION_H
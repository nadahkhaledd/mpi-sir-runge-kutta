#ifndef GRIDSIMULATION_H
#define GRIDSIMULATION_H

#include <vector>
#include <unordered_map>
#include <map>
#include <list>
#include <string>
#include "SIRCell.h"
#include "SIRModel.h"
#include <unordered_map>

class GridSimulation {
private:
    std::vector<SIRCell> grid;
    SIRModel model;
    int rank, size;
    std::unordered_map<int, std::vector<int>> neighborMap;
    std::unordered_map<int, int> cellIdToLocalIndex;
    std::unordered_map<int, std::vector<int>> ghostNeighborMap;

    
public:
    GridSimulation(const SIRModel& m, int mpiRank, int mpiSize);
    
    void setGrid(const std::vector<SIRCell>& initialGrid);
    
    std::vector<SIRCell>& getGrid();
    
    void exchangeGhostCells();
    
    int getLocalSize() const;
    
    void updateGrid();
    void updateGridNew();
    void setNeighborMap(const std::unordered_map<int, std::vector<int>>& map);
    void setGhostNeighborMap(const std::unordered_map<int, std::vector<int>>& map);

    
    std::vector<std::vector<double>> runSimulation();

    static std::map<std::string, int> createCellsMap();
    static std::map<int, std::list<int>> divideIntoBlocks(
        const std::map<std::string, int>& cells, int blockSize);
};

#endif // GRIDSIMULATION_H
#include <mpi.h>
#include "header/MPIHandler.h"
#include "header/SIRModel.h"
#include "header/CSVParser.h"
#include "header/GridSimulation.h"
#include <iostream>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

std::unordered_map<int, std::vector<int>> build2DGridNeighborMap(
    int rows, int cols,
    const std::unordered_map<int, int> &cellToBlock,
    std::unordered_map<int, std::vector<int>> &ghostNeighbors)
{
    std::unordered_map<int, std::vector<int>> neighbors;

    int totalCells = rows * cols;
    for (int i = 0; i < totalCells; ++i)
    {
        if (cellToBlock.find(i) == cellToBlock.end())
            continue;

        int row = i / cols;
        int col = i % cols;

        std::vector<int> neighborList;
        std::vector<int> directions = {
            (row > 0) ? (i - cols) : -1,
            (row < rows - 1) ? (i + cols) : -1,
            (col > 0) ? (i - 1) : -1,
            (col < cols - 1) ? (i + 1) : -1};

        for (int neighborId : directions)
        {
            if (neighborId != -1 && cellToBlock.find(neighborId) != cellToBlock.end())
            {
                if (cellToBlock.at(i) == cellToBlock.at(neighborId))
                {
                    neighborList.push_back(neighborId); // local
                }
                else
                {
                    ghostNeighbors[i].push_back(neighborId); // ghost
                    neighborList.push_back(neighborId);      // <-- ALSO keep it in the main list
                }
            }
        }

        neighbors[i] = neighborList;
    }

    return neighbors;
}

std::unordered_map<int, std::vector<int>> buildBlockNeighborMap(
    const std::map<int, std::list<int>> &allBlocks,
    const std::unordered_map<int, std::vector<int>> &cellNeighborMap)
{

    std::unordered_map<int, std::vector<int>> blockNeighborMap;
    if (allBlocks.empty() || cellNeighborMap.empty())
    {
        return blockNeighborMap;
    }

    std::unordered_map<int, int> cellToBlockMap;
    for (const auto &[blockId, cellList] : allBlocks)
    {
        for (int cellId : cellList)
        {
            cellToBlockMap[cellId] = blockId;
        }
    }

    for (const auto &[blockId, cellList] : allBlocks)
    {
        std::unordered_set<int> neighborBlockIdsSet;
        for (int cellId : cellList)
        {
            auto it_neighbors = cellNeighborMap.find(cellId);
            if (it_neighbors != cellNeighborMap.end())
            {
                for (int neighborCellId : it_neighbors->second)
                {
                    auto it_block = cellToBlockMap.find(neighborCellId);
                    if (it_block != cellToBlockMap.end())
                    {
                        int neighborBlockId = it_block->second;
                        if (neighborBlockId != blockId)
                        {
                            neighborBlockIdsSet.insert(neighborBlockId);
                        }
                    }
                }
            }
        }
        blockNeighborMap[blockId] = std::vector<int>(neighborBlockIdsSet.begin(), neighborBlockIdsSet.end());
        std::cerr << "DEBUG: Block " << blockId << " has " << blockNeighborMap[blockId].size() << " neighbors." << std::endl;
    }

    std::cerr << "DEBUG: Finished buildBlockNeighborMap. Size = " << blockNeighborMap.size() << std::endl;
    return blockNeighborMap;
}

int main(int argc, char *argv[])
{
    MPIHandler mpi(argc, argv);

    std::vector<std::vector<double>> fullData;
    std::map<int, std::list<int>> allBlocks;
    std::map<std::string, int> cells;
    std::unordered_map<int, std::vector<int>> blockNeighborMap;
    std::unordered_map<int, int> cellToBlock;

    int rows = 0, cols = 0;
    std::unordered_map<int, std::vector<int>> cellNeighborMap;
    std::unordered_map<int, std::vector<int>> ghostNeighborMap;

    if (mpi.getRank() == 0)
    {
        fullData = CSVParser::loadUSStateData("./data/sorted_initial_conditions.csv");
        if (fullData.empty())
        {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        cells = GridSimulation::createCellsMap();
        allBlocks = GridSimulation::divideIntoOptimalBlocks(cells, mpi.getSize());

        for (const auto &[blockId, cellList] : allBlocks)
        {
            for (int cell : cellList)
            {
                cellToBlock[cell] = blockId;
            }
        }

        int totalCells = fullData.size();
        cols = std::ceil(std::sqrt(totalCells));
        rows = (totalCells + cols - 1) / cols;
        std::cout << "Adjusted grid: rows = " << rows << ", cols = " << cols << "\n";

        cellNeighborMap = build2DGridNeighborMap(rows, cols, cellToBlock, ghostNeighborMap);

        std::cout << "Rank " << mpi.getRank() << " Ghost Neighbors (if any):\n";
        for (const auto &[cellId, ghosts] : ghostNeighborMap)
        {
            std::cout << "  Cell " << cellId << " has ghost neighbors: ";
            for (int ghost : ghosts)
                std::cout << ghost << " ";
            std::cout << "\n";
        }

        blockNeighborMap = buildBlockNeighborMap(allBlocks, cellNeighborMap);
        std::cout << "DEBUG: blockNeighborMap size = " << blockNeighborMap.size() << "\n";

        for (const auto &[blockId, neighborBlocks] : blockNeighborMap)
        {
            std::cout << "Block " << blockId << " neighbors: ";
            for (int neighborBlock : neighborBlocks)
            {
                std::cout << neighborBlock << " ";
            }
            std::cout << "\n";
        }
    }

    std::map<int, std::list<int>> localBlocks = mpi.distributeBlocks(allBlocks);
    std::map<int, std::vector<double>> localCellData = mpi.getDataForLocalBlocks(localBlocks, fullData);
    blockNeighborMap = mpi.broadcastBlockNeighborMap(blockNeighborMap);

    SIRModel model(0.3, 0.1, 0.2, 100);
    GridSimulation simulation(model, mpi.getRank(), mpi.getSize());

    simulation.setCellNeighborMap(cellNeighborMap);
    simulation.setGhostNeighborMap(ghostNeighborMap);

    std::unordered_map<int, int> blockToRankMap;
    if (mpi.getRank() == 0)
    {
        int totalBlocks = allBlocks.size();
        int blockCounter = 0;
        for (const auto &[blockId, _] : allBlocks)
        {
            int targetRank = -1;
            int tempBlocksPerProc = totalBlocks / mpi.getSize();
            int tempExtraBlocks = totalBlocks % mpi.getSize();
            for (int r = 0; r < mpi.getSize(); ++r)
            {
                int rNumBlocks = (r < tempExtraBlocks) ? (tempBlocksPerProc + 1) : tempBlocksPerProc;
                int rankStartBlockIndex = (r < tempExtraBlocks) ? (r * (tempBlocksPerProc + 1)) : (r * tempBlocksPerProc + tempExtraBlocks);
                if (blockCounter >= rankStartBlockIndex && blockCounter < rankStartBlockIndex + rNumBlocks)
                {
                    targetRank = r;
                    break;
                }
            }
            if (targetRank != -1)
            {
                blockToRankMap[blockId] = targetRank;
            }
            blockCounter++;
        }
    }

    int map_size = (mpi.getRank() == 0) ? blockToRankMap.size() : 0;
    MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> map_data(map_size * 2);
    if (mpi.getRank() == 0)
    {
        int i = 0;
        for (const auto &pair : blockToRankMap)
        {
            map_data[i++] = pair.first;
            map_data[i++] = pair.second;
        }
    }
    MPI_Bcast(map_data.data(), map_size * 2, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi.getRank() != 0)
    {
        for (int i = 0; i < map_size * 2; i += 2)
        {
            blockToRankMap[map_data[i]] = map_data[i + 1];
        }
    }
    simulation.setBlockToRankMap(blockToRankMap);
    simulation.setGridFromLocalData(localBlocks, localCellData);
    simulation.setBlockInfo(localBlocks, blockNeighborMap);

    std::vector<std::vector<double>> localResults = simulation.runSimulation();
    std::vector<int> recvCounts;      // To hold the number of doubles each rank contributes
    std::vector<int> displacements;   // To hold the offset into the gathered buffer per rank

    std::vector<double> globalResults = mpi.gatherResults(localResults, recvCounts, displacements);

    mpi.writeResults(globalResults, recvCounts, displacements);

    return 0;
}

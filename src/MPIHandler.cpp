#include "../header/MPIHandler.h"
#include "../header/CSVParser.h" // Needed for mapToSIR
#include "../header/SIRCell.h"   // Needed for SIRCell type

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>
#include <unordered_map>
#include <cstring>
#include <stdexcept>
#include <set>

// Constructor: Initialize MPI
MPIHandler::MPIHandler(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

// Destructor: Finalize MPI
MPIHandler::~MPIHandler() {
    int finalized_flag;
    MPI_Finalized(&finalized_flag);
    if (!finalized_flag) {
        MPI_Finalize();
    }
}

// Getters
int MPIHandler::getRank() const {
    return rank;
}

int MPIHandler::getSize() const {
    return size;
}

// gatherResults: Gather [time, avgS, avgI, avgR] from each process
std::vector<double> MPIHandler::gatherResults(const std::vector<std::vector<double>>& localResults,
                                              std::vector<int>& outCounts,
                                              std::vector<int>& outDispls) {
    int doublesPerStep = 4;
    int localSteps = localResults.size();
    int localDataSize = localSteps * doublesPerStep;

    std::vector<double> localFlat;
    localFlat.reserve(localDataSize);
    for (const auto& stepData : localResults) {
        if (stepData.size() == static_cast<size_t>(doublesPerStep)) {
            localFlat.insert(localFlat.end(), stepData.begin(), stepData.end());
        } else {
            localFlat.insert(localFlat.end(), doublesPerStep, 0.0);
        }
    }
    localDataSize = localFlat.size();

    std::vector<int> recvCounts(size);
    std::vector<int> displacements(size);
    MPI_Gather(&localDataSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<double> globalFlat;
    if (rank == 0) {
        displacements[0] = 0;
        for (int i = 1; i < size; ++i)
            displacements[i] = displacements[i - 1] + recvCounts[i - 1];
        globalFlat.resize(displacements[size - 1] + recvCounts[size - 1]);
    }

    MPI_Gatherv(localFlat.data(), localDataSize, MPI_DOUBLE,
                globalFlat.data(), recvCounts.data(), displacements.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    outCounts = recvCounts;
    outDispls = displacements;

    if (rank == 0) {
        for (size_t i = 0; i < globalFlat.size(); i += 4) {
            double S = globalFlat[i + 1];
            double I = globalFlat[i + 2];
            double R = globalFlat[i + 3];
            double sum = S + I + R;

            if (sum > 0) {
                globalFlat[i + 1] = std::max(0.0, std::min(1.0, S / sum));
                globalFlat[i + 2] = std::max(0.0, std::min(1.0, I / sum));
                globalFlat[i + 3] = std::max(0.0, std::min(1.0, R / sum));
            } else {
                globalFlat[i + 1] = 1.0;
                globalFlat[i + 2] = 0.0;
                globalFlat[i + 3] = 0.0;
            }
        }
    }
    return globalFlat;
}

// writeResults: now receives correct counts/displacements
void MPIHandler::writeResults(const std::vector<double>& globalFlat,
                              const std::vector<int>& recvCounts,
                              const std::vector<int>& displacements) {
    if (rank == 0) {
        std::ofstream outfile("./data/simulation_results.csv");
        if (!outfile) {
            return;
        }

        outfile << "Rank,Time,S_avg,I_avg,R_avg\n";
        int doublesPerStep = 4;
        for (int proc = 0; proc < size; ++proc) {
            int startIdx = displacements[proc];
            int numDoubles = recvCounts[proc];
            if (numDoubles % doublesPerStep != 0) {
                continue;
            }
            int numSteps = numDoubles / doublesPerStep;
            for (int i = 0; i < numSteps; ++i) {
                int idx = startIdx + i * doublesPerStep;
                outfile << proc << ","
                        << globalFlat[idx + 0] << ","
                        << globalFlat[idx + 1] << ","
                        << globalFlat[idx + 2] << ","
                        << globalFlat[idx + 3] << "\n";
            }
        }
        outfile.close();
    }
}

// distributeBlocks: Distributes ONLY block structure (IDs)
std::map<int, std::list<int>> MPIHandler::distributeBlocks(
    const std::map<int, std::list<int>>& allBlocks) {
    std::cout << "Rank " << rank << ": Starting block distribution...\n";

    std::map<int, std::list<int>> localBlocks;
    int totalBlocks = 0;

    if (rank == 0) {
        totalBlocks = allBlocks.size();
        std::cout << "Rank 0: Distributing " << totalBlocks << " blocks among " << size << " processes.\n";
    }
    MPI_Bcast(&totalBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (totalBlocks == 0) {
        return localBlocks;
    }

    int blocksPerProc = totalBlocks / size;
    int extraBlocks = totalBlocks % size;
    int numMyBlocks = (rank < extraBlocks) ? (blocksPerProc + 1) : blocksPerProc;
    int startBlockIndex = (rank < extraBlocks) ? (rank * (blocksPerProc + 1)) : (rank * blocksPerProc + extraBlocks);

    if (rank == 0) {
        int currentBlockIndex = 0;
        for (const auto& [blockId, cellList] : allBlocks) {
            if (currentBlockIndex >= startBlockIndex && currentBlockIndex < startBlockIndex + numMyBlocks) {
                localBlocks[blockId] = cellList;
            }
            currentBlockIndex++;
        }

        int blockCounter = 0;
        for (const auto& [blockId, cellList] : allBlocks) {
            int targetRank = -1;
            int tempBlocksPerProc = totalBlocks / size;
            int tempExtraBlocks = totalBlocks % size;
            for (int r = 0; r < size; ++r) {
                int rNumBlocks = (r < tempExtraBlocks) ? (tempBlocksPerProc + 1) : tempBlocksPerProc;
                int rankStartBlockIndex = (r < tempExtraBlocks) ? (r * (tempBlocksPerProc + 1)) : (r * tempBlocksPerProc + tempExtraBlocks);
                if (blockCounter >= rankStartBlockIndex && blockCounter < rankStartBlockIndex + rNumBlocks) {
                    targetRank = r;
                    break;
                }
            }

            if (targetRank > 0 && targetRank < size) {
                std::vector<int> blockData;
                blockData.push_back(blockId);
                blockData.push_back(static_cast<int>(cellList.size()));
                blockData.insert(blockData.end(), cellList.begin(), cellList.end());

                int dataSize = blockData.size();
                MPI_Send(&dataSize, 1, MPI_INT, targetRank, 0, MPI_COMM_WORLD);
                if (dataSize > 0) {
                    MPI_Send(blockData.data(), dataSize, MPI_INT, targetRank, 1, MPI_COMM_WORLD);
                }
            }
            blockCounter++;
        }

        int terminateSignal = -1;
        for (int proc = 1; proc < size; ++proc) {
            MPI_Send(&terminateSignal, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }

    } else {
        while (true) {
            int dataSize;
            MPI_Status status;
            MPI_Recv(&dataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

            if (dataSize == -1) {
                break;
            }

            if (dataSize > 0) {
                std::vector<int> blockData(dataSize);
                MPI_Recv(blockData.data(), dataSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (blockData.size() >= 2) {
                    int blockId = blockData[0];
                    int numCells = blockData[1];
                    if (blockData.size() == static_cast<size_t>(2 + numCells)) {
                        std::list<int> cellList;
                        for (int i = 0; i < numCells; ++i) {
                            cellList.push_back(blockData[2 + i]);
                        }
                        localBlocks[blockId] = cellList;
                    }
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Rank " << rank << ": Finished block distribution. Local blocks: " << localBlocks.size() << "\n";
    return localBlocks;
}

// getDataForLocalBlocks: Fetches the necessary initial condition data rows
std::map<int, std::vector<double>> MPIHandler::getDataForLocalBlocks(
    const std::map<int, std::list<int>>& localBlocks,
    const std::vector<std::vector<double>>& fullData) {
    std::cout << "Rank " << rank << ": Starting data distribution for local blocks...\n";

    std::map<int, std::vector<double>> localCellData;

    std::set<int> neededCellIdsSet;
    for (const auto& [blockId, cellList] : localBlocks) {
        neededCellIdsSet.insert(cellList.begin(), cellList.end());
    }
    std::vector<int> neededCellIds(neededCellIdsSet.begin(), neededCellIdsSet.end());

    if (rank == 0) {
        int myRequestSize = neededCellIds.size();
        std::vector<int> requestSizes(size);
        MPI_Gather(&myRequestSize, 1, MPI_INT, requestSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<int> requestDispls(size);
        int totalRequestedIds = 0;
        requestDispls[0] = 0;
        totalRequestedIds = requestSizes[0];
        for (int i = 1; i < size; ++i) {
            requestDispls[i] = requestDispls[i - 1] + requestSizes[i - 1];
            totalRequestedIds += requestSizes[i];
        }

        std::vector<int> gatheredIdsBuffer(totalRequestedIds);
        MPI_Gatherv(neededCellIds.data(), myRequestSize, MPI_INT,
                    gatheredIdsBuffer.data(), requestSizes.data(), requestDispls.data(), MPI_INT,
                    0, MPI_COMM_WORLD);

        int doublesPerCell = 0;
        if (!fullData.empty() && !fullData[0].empty()) {
            doublesPerCell = fullData[0].size();
        } else if (totalRequestedIds > 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::vector<int> sendDataSizes(size);
        std::vector<int> sendDataDispls(size);
        std::vector<double> flatSendDataBuffer;

        int currentSendDispl = 0;
        for (int targetRank = 0; targetRank < size; ++targetRank) {
            int startIdx = requestDispls[targetRank];
            int numIds = requestSizes[targetRank];
            std::vector<double> rankDataBuffer;

            for (int k = 0; k < numIds; ++k) {
                int cellId = gatheredIdsBuffer[startIdx + k];
                if (cellId >= 0 && static_cast<size_t>(cellId) < fullData.size()) {
                    const auto& rowData = fullData[cellId];
                    rankDataBuffer.insert(rankDataBuffer.end(), rowData.begin(), rowData.end());

                    if (targetRank == 0) {
                        localCellData[cellId] = rowData;
                    }
                }
            }
            sendDataSizes[targetRank] = rankDataBuffer.size();
            sendDataDispls[targetRank] = currentSendDispl;
            flatSendDataBuffer.insert(flatSendDataBuffer.end(), rankDataBuffer.begin(), rankDataBuffer.end());
            currentSendDispl += rankDataBuffer.size();
        }

        MPI_Bcast(&doublesPerCell, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(sendDataSizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

        std::vector<double> localRecvBuffer_Rank0(sendDataSizes[0]);
        MPI_Scatterv(flatSendDataBuffer.data(), sendDataSizes.data(), sendDataDispls.data(), MPI_DOUBLE,
                     localRecvBuffer_Rank0.data(), sendDataSizes[0], MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

    } else {
        int myRequestSize = neededCellIds.size();
        MPI_Gather(&myRequestSize, 1, MPI_INT, nullptr, 0, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Gatherv(neededCellIds.data(), myRequestSize, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT,
                    0, MPI_COMM_WORLD);

        int doublesPerCell = 0;
        MPI_Bcast(&doublesPerCell, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (doublesPerCell <= 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::vector<int> dataRecvSizes(size);
        MPI_Bcast(dataRecvSizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

        int myRecvSize = dataRecvSizes[rank];
        std::vector<double> localRecvBuffer(myRecvSize);

        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE,
                     localRecvBuffer.data(), myRecvSize, MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        if (myRecvSize % doublesPerCell != 0 && myRecvSize != 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            size_t numCellsReceived = (myRecvSize == 0) ? 0 : static_cast<size_t>(myRecvSize / doublesPerCell);
            size_t bufferIdx = 0;

            for (size_t i = 0; i < numCellsReceived; ++i) {
                int cellId = neededCellIds[i];
                std::vector<double> rowData(localRecvBuffer.begin() + bufferIdx,
                                            localRecvBuffer.begin() + bufferIdx + doublesPerCell);
                localCellData[cellId] = rowData;
                bufferIdx += doublesPerCell;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Rank " << rank << ": Received data for " << localCellData.size() << " cells.\n";
    return localCellData;
}

// broadcastBlockNeighborMap: Broadcasts the block adjacency map
std::unordered_map<int, std::vector<int>> MPIHandler::broadcastBlockNeighborMap(
    const std::unordered_map<int, std::vector<int>>& mapToSend) {
    std::cout << "Rank " << rank << ": Broadcasting block neighbor map...\n";

    long long totalBytes = 0;
    std::vector<char> buffer;
    if (rank == 0) {
        std::vector<int> flatIntBuffer;
        flatIntBuffer.push_back(static_cast<int>(mapToSend.size()));
        for (const auto& [key, vec] : mapToSend) {
            flatIntBuffer.push_back(key);
            flatIntBuffer.push_back(static_cast<int>(vec.size()));
            flatIntBuffer.insert(flatIntBuffer.end(), vec.begin(), vec.end());
        }
        totalBytes = static_cast<long long>(flatIntBuffer.size()) * sizeof(int);
        buffer.resize(totalBytes);
        if (totalBytes > 0) {
            memcpy(buffer.data(), flatIntBuffer.data(), totalBytes);
        }
    }
    MPI_Bcast(&totalBytes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    if (totalBytes == 0) {
        return (rank == 0) ? mapToSend : std::unordered_map<int, std::vector<int>>{};
    }
    if (rank != 0) {
        buffer.resize(totalBytes);
    }
    if (totalBytes > 0) {
        MPI_Bcast(buffer.data(), totalBytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    std::unordered_map<int, std::vector<int>> receivedMap;
    if (rank != 0) {
        size_t numInts = static_cast<size_t>(totalBytes) / sizeof(int);
        std::vector<int> flatIntBuffer(numInts);
        if (totalBytes > 0) memcpy(flatIntBuffer.data(), buffer.data(), totalBytes);
        int numEntries = flatIntBuffer[0];
        size_t currentIndex = 1;
        receivedMap.reserve(numEntries);
        for (int i = 0; i < numEntries; ++i) {
            int key = flatIntBuffer[currentIndex++];
            int numNeighbors = flatIntBuffer[currentIndex++];
            std::vector<int> neighbors;
            if (numNeighbors > 0) {
                neighbors.assign(flatIntBuffer.begin() + currentIndex, flatIntBuffer.begin() + currentIndex + numNeighbors);
                currentIndex += numNeighbors;
            }
            receivedMap[key] = std::move(neighbors);
        }
    }
    std::cout << "Rank " << rank << ": Received block neighbor map with " << receivedMap.size() << " entries.\n";
    return (rank == 0) ? mapToSend : receivedMap;
}
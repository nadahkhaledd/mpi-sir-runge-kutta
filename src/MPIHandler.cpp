#include "../header/MPIHandler.h"
#include "../header/CSVParser.h" // Needed for mapToSIR
#include "../header/SIRCell.h"   // Needed for SIRCell type
#include "../header/TimingUtils.h"
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
#include <iomanip>      // For std::fixed, std::setprecision
#include <string>       // For std::string
#include <numeric>      // For std::accumulate (not used for summary here)
#include <algorithm>    // For std::max/min in gatherResults, and MPI_MAX reduce op

// Constructor: Initialize MPI
MPIHandler::MPIHandler(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // std::cout << "Rank " << rank << ": MPI initialized by MPIHandler." << std::endl;
}

// Destructor: Finalize MPI
MPIHandler::~MPIHandler() {
    // std::cout << "Rank " << rank << ": MPIHandler finalizing MPI." << std::endl;
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
                                              std::vector<int>& outCounts, // Populated on rank 0
                                              std::vector<int>& outDispls) { // Populated on rank 0
    // std::cout << "Rank " << rank << ": MPIHandler::gatherResults called." << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    double totalPhaseStartTime = MPI_Wtime();

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

    if (rank == 0) {
        outCounts.assign(size, 0);
        outDispls.assign(size, 0);
    }
    std::vector<int> recvCounts(size);
    std::vector<int> displacements(size);

    MPI_Barrier(MPI_COMM_WORLD);
    double gatherSizesStartTime = MPI_Wtime();
    MPI_Gather(&localDataSize, 1, MPI_INT, recvCounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6) << "TIMING [Rank 0, Phase: gatherResults_GatherSizes]: " << (MPI_Wtime() - gatherSizesStartTime) << "s\n";
    }

    std::vector<double> globalFlat;
    int totalGlobalSize = 0;
    if (rank == 0) {
        displacements[0] = 0;
        totalGlobalSize = recvCounts[0];
        for (int i = 1; i < size; ++i) {
            if (i < static_cast<int>(recvCounts.size()) && (i - 1) < static_cast<int>(recvCounts.size())) {
                 displacements[i] = displacements[i - 1] + recvCounts[i - 1];
                 totalGlobalSize += recvCounts[i];
            }
        }
        if (totalGlobalSize > 0) globalFlat.resize(totalGlobalSize);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double gatherVDataStartTime = MPI_Wtime();
    MPI_Gatherv(localFlat.data(), localDataSize, MPI_DOUBLE,
                globalFlat.data(), recvCounts.data(), displacements.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6) << "TIMING [Rank 0, Phase: gatherResults_GatherVData]: " << (MPI_Wtime() - gatherVDataStartTime) << "s\n";
    }

    if (rank == 0) {
        outCounts = recvCounts;
        outDispls = displacements;
    }

    if (rank == 0) {
        double normStartTime = MPI_Wtime();
        for (size_t i = 0; i < globalFlat.size(); i += 4) {
            if (i + 3 >= globalFlat.size()) break;
            double S = globalFlat[i + 1]; double I = globalFlat[i + 2]; double R = globalFlat[i + 3];
            double sum = S + I + R;
            if (sum > 1e-9) {
                globalFlat[i + 1] = std::max(0.0, std::min(1.0, S / sum));
                globalFlat[i + 2] = std::max(0.0, std::min(1.0, I / sum));
                globalFlat[i + 3] = std::max(0.0, std::min(1.0, R / sum));
            } else { globalFlat[i + 1] = 1.0; globalFlat[i + 2] = 0.0; globalFlat[i + 3] = 0.0; }
        }
        std::cout << std::fixed << std::setprecision(6) << "TIMING [Rank 0, Phase: gatherResults_Normalization]: " << (MPI_Wtime() - normStartTime) << "s\n";
    }

    double phaseDuration = MPI_Wtime() - totalPhaseStartTime;
    double maxDuration;
    MPI_Reduce(&phaseDuration, &maxDuration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 
    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: gatherResults_Total_Max]: " << maxDuration << "s\n";
    }
    return globalFlat;
}

// writeResults: now receives correct counts/displacements
void MPIHandler::writeResults(const std::vector<double>& globalFlat,
                              const std::vector<int>& recvCounts,
                              const std::vector<int>& displacements) {
    double ioStartTime = 0.0;
    if (rank == 0) {
        ioStartTime = MPI_Wtime();
        if (!globalFlat.empty() && (recvCounts.empty() || displacements.empty() || recvCounts.size()!=static_cast<size_t>(size) || displacements.size()!=static_cast<size_t>(size))) {
             std::cerr << "Rank 0 Error in writeResults: recvCounts/displacements invalid." << std::endl; return;
        }
        std::ofstream outfile("./data/simulation_results.csv");
        if (!outfile) { std::cerr << "Rank 0 Error: Cannot open ./data/simulation_results.csv" << std::endl; return; }
        outfile << "Rank,Time,S_avg,I_avg,R_avg\n";
        int doublesPerStep = 4;
        for (int proc = 0; proc < size; ++proc) {
            if (proc >= static_cast<int>(displacements.size()) || proc >= static_cast<int>(recvCounts.size())) continue;
            int startIdx = displacements[proc]; int numDoubles = recvCounts[proc];
            if (numDoubles < 0 || (numDoubles % doublesPerStep != 0 && numDoubles !=0) ) continue;
            if (numDoubles == 0) continue;
            int numSteps = numDoubles / doublesPerStep;
            for (int i = 0; i < numSteps; ++i) {
                int idx = startIdx + i * doublesPerStep;
                if (idx + doublesPerStep - 1 < static_cast<int>(globalFlat.size())) {
                    outfile << proc << "," << globalFlat[idx + 0] << "," << globalFlat[idx + 1] << "," << globalFlat[idx + 2] << "," << globalFlat[idx + 3] << "\n";
                } else { break; }
            }
        }
        outfile.close();
        std::cout << "Rank 0: Results written to ./data/simulation_results.csv" << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << std::fixed << std::setprecision(6) << "TIMING [Rank 0, Phase: writeResults_FileIO]: " << (MPI_Wtime() - ioStartTime) << "s\n";
    }
}
// distributeBlocks: Distributes ONLY block structure (IDs)
std::map<int, std::list<int>> MPIHandler::distributeBlocks(
    const std::map<int, std::list<int>>& allBlocks) {

    // General informational message (optional, printed by each rank before timing)
    // std::cout << "Rank " << rank << ": MPIHandler::distributeBlocks called." << std::endl;

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting overall phase timing
    double totalPhase_StartTime = TimingUtils::startTimer(); // Start timer for the whole function

    std::map<int, std::list<int>> localBlocks;
    int totalBlocks = 0;

    // --- Timing MPI_Bcast of totalBlocks ---
    MPI_Barrier(MPI_COMM_WORLD); // Sync for this specific collective's timing start
    double bcastTotalBlocks_StartTime = TimingUtils::startTimer();
    if (rank == 0) {
        totalBlocks = allBlocks.size();
        // Informational print by Rank 0
        std::cout << "Rank 0 (distributeBlocks): Distributing " << totalBlocks << " blocks among " << size << " processes.\n";
    }
    MPI_Bcast(&totalBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // Ensure Bcast completes on all before Rank 0 measures/prints
    // Rank 0 prints its observed time for this collective, and logs it.
    TimingUtils::printAndLogRank0PhaseTime("distributeBlocks_Bcast_TotalBlocks", bcastTotalBlocks_StartTime, rank);
    // --- End Timing Bcast ---


    if (totalBlocks == 0) {
        // Optional informational print
        // if (rank == 0) std::cout << "Rank 0 (distributeBlocks): No blocks to distribute." << std::endl;
        double totalEmpty_LocalDuration = TimingUtils::stopTimer(totalPhase_StartTime);
        // All ranks participate, Rank 0 prints summary and logs it.
        TimingUtils::printAndLogTimingSummary("distributeBlocks_Total_Empty", totalEmpty_LocalDuration, rank, size);
        return localBlocks;
    }

    // --- Block assignment logic (unchanged from your working version) ---
    int blocksPerProc = totalBlocks / size;
    int extraBlocks = totalBlocks % size;
    int numMyBlocks = (rank < extraBlocks) ? (blocksPerProc + 1) : blocksPerProc;
    int startBlockIndex = (rank < extraBlocks) ? (rank * (blocksPerProc + 1)) : (rank * blocksPerProc + extraBlocks);

    // --- Timing Communication Loop (Send/Recv of block structures) ---
    MPI_Barrier(MPI_COMM_WORLD); // Sync before the parallel send/recv loop
    double commLoop_StartTime = TimingUtils::startTimer(); // Each rank starts its own timer for its part
    if (rank == 0) {
        int currentBlockIndex = 0;
        for (const auto& p : allBlocks) { // Using p for std::pair
            if (currentBlockIndex >= startBlockIndex && currentBlockIndex < startBlockIndex + numMyBlocks) {
                localBlocks[p.first] = p.second;
            }
            currentBlockIndex++;
        }
        int blockCounter = 0;
        for (const auto& p : allBlocks) { // Using p for std::pair
            int targetRank = -1;
            // Determine targetRank based on blockCounter (your existing logic)
            int tempBlocksPerProc = totalBlocks / size; // Recalculate for safety or pass from above
            int tempExtraBlocks = totalBlocks % size;
            for (int r = 0; r < size; ++r) {
                int rNB = (r < tempExtraBlocks) ? (tempBlocksPerProc + 1) : tempBlocksPerProc;
                int rS = (r < tempExtraBlocks) ? (r * (tempBlocksPerProc + 1)) : (r * tempBlocksPerProc + tempExtraBlocks);
                if (blockCounter >= rS && blockCounter < rS + rNB) {
                    targetRank = r;
                    break;
                }
            }
            if (targetRank > 0 && targetRank < size) { // Send only to other valid ranks
                std::vector<int> bd; // blockData
                bd.push_back(p.first); // blockId
                bd.push_back(static_cast<int>(p.second.size())); // numCells
                bd.insert(bd.end(), p.second.begin(), p.second.end()); // cellList
                int ds = bd.size(); // dataSize
                MPI_Send(&ds, 1, MPI_INT, targetRank, 0, MPI_COMM_WORLD);
                if (ds > 0) MPI_Send(bd.data(), ds, MPI_INT, targetRank, 1, MPI_COMM_WORLD);
            }
            blockCounter++;
        }
        int termSig = -1; // terminateSignal
        for (int proc_idx = 1; proc_idx < size; ++proc_idx) MPI_Send(&termSig, 1, MPI_INT, proc_idx, 0, MPI_COMM_WORLD);
    } else { // Other ranks receive
        while (true) {
            int ds; // dataSize
            MPI_Status st;
            MPI_Recv(&ds, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
            if (ds == -1) break; // Termination signal
            if (ds > 0) {
                std::vector<int> bd(ds); // blockData
                MPI_Recv(bd.data(), ds, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (bd.size() >= 2) { // At least blockId and numCells
                    int bId = bd[0]; // blockId
                    int nC = bd[1];  // numCells
                    if (bd.size() == static_cast<size_t>(2 + nC)) {
                        std::list<int> cl; // cellList
                        for (int i = 0; i < nC; ++i) cl.push_back(bd[2 + i]);
                        localBlocks[bId] = cl;
                    } else {
                        // Error: size mismatch
                        std::cerr << "Rank " << rank << " Error in distributeBlocks: Received block data size mismatch for block " << bId << std::endl;
                    }
                } else {
                    // Error: received data too small
                     std::cerr << "Rank " << rank << " Error in distributeBlocks: Received block data too small (size " << bd.size() << ")" << std::endl;
                }
            } else if (ds < -1 ) { // ds == 0 might be possible if a block has 0 cells, but less likely
                 std::cerr << "Rank " << rank << " Warning/Error in distributeBlocks: Received unusual dataSize " << ds << std::endl;
            }
        }
    }
    double commLoop_LocalDuration = TimingUtils::stopTimer(commLoop_StartTime); // Each rank stops its timer
    // Gather all local durations for the CommLoop and print Min/Max/Avg by Rank 0, also logs it
    TimingUtils::printAndLogTimingSummary("distributeBlocks_CommLoop", commLoop_LocalDuration, rank, size);
    // --- End Timing Communication Loop ---


    double total_LocalDuration = TimingUtils::stopTimer(totalPhase_StartTime); // Each rank gets its total duration for the function
    // Gather all total durations and print Min/Max/Avg by Rank 0, also logs it
    TimingUtils::printAndLogTimingSummary("distributeBlocks_Total", total_LocalDuration, rank, size);

    // General informational message (optional, printed by each rank after timing)
    // std::cout << "Rank " << rank << ": Finished block distribution. Local blocks: " << localBlocks.size() << "\n";
    return localBlocks;
}

// getDataForLocalBlocks: Fetches the necessary initial condition data rows
std::map<int, std::vector<double>> MPIHandler::getDataForLocalBlocks(
    const std::map<int, std::list<int>>& localBlocks,
    const std::vector<std::vector<double>>& fullData) {
    // std::cout << "Rank " << rank << ": Starting data distribution for local blocks...\n"; // Moved
    MPI_Barrier(MPI_COMM_WORLD);
    double totalPhaseStartTime = MPI_Wtime();

    std::map<int, std::vector<double>> localCellData;
    std::set<int> neededCellIdsSet;
    for (const auto& [blockId, cellList] : localBlocks) {
        neededCellIdsSet.insert(cellList.begin(), cellList.end());
    }
    std::vector<int> neededCellIds(neededCellIdsSet.begin(), neededCellIdsSet.end());

    MPI_Barrier(MPI_COMM_WORLD);
    double gatherRequestsStartTime = MPI_Wtime();
    std::vector<int> requestSizes(size);
    std::vector<int> requestDispls(size);
    std::vector<int> gatheredIdsBuffer; // Will be resized on rank 0

    int myRequestSize = neededCellIds.size();
    MPI_Gather(&myRequestSize, 1, MPI_INT, requestSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        int totalRequestedIds = 0; requestDispls[0] = 0; totalRequestedIds = requestSizes[0];
        for (int i = 1; i < size; ++i) { requestDispls[i] = requestDispls[i-1] + requestSizes[i-1]; totalRequestedIds += requestSizes[i]; }
        if(totalRequestedIds > 0) gatheredIdsBuffer.resize(totalRequestedIds);
    }
    MPI_Gatherv(neededCellIds.data(), myRequestSize, MPI_INT,
                gatheredIdsBuffer.data(), requestSizes.data(), requestDispls.data(), MPI_INT,
                0, MPI_COMM_WORLD);
    double gatherReqDuration = MPI_Wtime() - gatherRequestsStartTime;
    double maxGatherReqDuration; MPI_Reduce(&gatherReqDuration, &maxGatherReqDuration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: getDataForLocalBlocks_GatherRequests_Max]: " << maxGatherReqDuration << "s\n";


    MPI_Barrier(MPI_COMM_WORLD);
    double dataDistCommStartTime = MPI_Wtime();
    int doublesPerCell = 0;
    std::vector<int> sendDataSizes(size); // For Scatterv counts

    if (rank == 0) {
        bool anyReq=false; for(int s:requestSizes)if(s>0)anyReq=true;
        if(!fullData.empty() && !fullData[0].empty())doublesPerCell=fullData[0].size();else if(anyReq)MPI_Abort(MPI_COMM_WORLD,1);
        std::vector<int>sdDispls(size); std::vector<double>flatSDBuf; int currSDispl=0;
        for(int tr=0;tr<size;++tr){int sIdx=requestDispls[tr];int nIds=requestSizes[tr];std::vector<double>rDBuf;for(int k=0;k<nIds;++k){int cId=gatheredIdsBuffer[sIdx+k];if(cId>=0&&static_cast<size_t>(cId)<fullData.size()){const auto&rd=fullData[cId];if(static_cast<int>(rd.size())!=doublesPerCell&&doublesPerCell>0)MPI_Abort(MPI_COMM_WORLD,1);rDBuf.insert(rDBuf.end(),rd.begin(),rd.end());if(tr==0)localCellData[cId]=rd;}} sendDataSizes[tr]=rDBuf.size();sdDispls[tr]=currSDispl;flatSDBuf.insert(flatSDBuf.end(),rDBuf.begin(),rDBuf.end());currSDispl+=rDBuf.size();}
        MPI_Bcast(&doublesPerCell, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(sendDataSizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
        std::vector<double>lrBufR0(sendDataSizes[0]>0?sendDataSizes[0]:0);MPI_Scatterv(flatSDBuf.data(),sendDataSizes.data(),sdDispls.data(),MPI_DOUBLE,(sendDataSizes[0]>0?lrBufR0.data():nullptr),sendDataSizes[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
    } else {
        MPI_Bcast(&doublesPerCell,1,MPI_INT,0,MPI_COMM_WORLD);if(doublesPerCell<=0&&neededCellIds.size()>0)MPI_Abort(MPI_COMM_WORLD,1);
        MPI_Bcast(sendDataSizes.data(),size,MPI_INT,0,MPI_COMM_WORLD);int myRS=(rank<static_cast<int>(sendDataSizes.size()))?sendDataSizes[rank]:0;std::vector<double>lrBuf(myRS);MPI_Scatterv(nullptr,nullptr,nullptr,MPI_DOUBLE,(myRS>0?lrBuf.data():nullptr),myRS,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(myRS>0){if(doublesPerCell == 0 && myRS > 0) MPI_Abort(MPI_COMM_WORLD,1); if(myRS % doublesPerCell != 0) MPI_Abort(MPI_COMM_WORLD,1);size_t nCRcvd=static_cast<size_t>(myRS/doublesPerCell);size_t bIdx=0;if(neededCellIds.size()<nCRcvd)MPI_Abort(MPI_COMM_WORLD,1);for(size_t i=0;i<nCRcvd;++i){if(bIdx+doublesPerCell>static_cast<size_t>(myRS))MPI_Abort(MPI_COMM_WORLD,1);int cId=neededCellIds[i];std::vector<double>rd(lrBuf.begin()+bIdx,lrBuf.begin()+bIdx+doublesPerCell);localCellData[cId]=rd;bIdx+=doublesPerCell;}}
    }
    double dataDistCommDuration = MPI_Wtime() - dataDistCommStartTime;
    double maxDataDistCommDuration; MPI_Reduce(&dataDistCommDuration, &maxDataDistCommDuration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: getDataForLocalBlocks_DataDistComm_Max]: " << maxDataDistCommDuration << "s\n";


    double totalDuration = MPI_Wtime() - totalPhaseStartTime;
    double maxTotalDurationData; MPI_Reduce(&totalDuration, &maxTotalDurationData, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: getDataForLocalBlocks_Total_Max]: " << maxTotalDurationData << "s\n";
    // std::cout << "Rank " << rank << ": Received data for " << localCellData.size() << " cells.\n"; // Moved
    return localCellData;
}

// broadcastBlockNeighborMap: Broadcasts the block adjacency map
std::unordered_map<int, std::vector<int>> MPIHandler::broadcastBlockNeighborMap(
    const std::unordered_map<int, std::vector<int>>& mapToSend) {
    // std::cout << "Rank " << rank << ": Broadcasting block neighbor map...\n"; // Moved
    MPI_Barrier(MPI_COMM_WORLD);
    double totalPhaseStartTime = MPI_Wtime();

    long long totalBytes = 0;
    std::vector<char> buffer;
    if (rank == 0) {
        std::vector<int> flatIntBuffer; flatIntBuffer.push_back(static_cast<int>(mapToSend.size()));
        for (const auto& [key, vec] : mapToSend) { flatIntBuffer.push_back(key); flatIntBuffer.push_back(static_cast<int>(vec.size())); flatIntBuffer.insert(flatIntBuffer.end(), vec.begin(), vec.end()); }
        totalBytes = static_cast<long long>(flatIntBuffer.size()) * sizeof(int);
        if(totalBytes > 0) buffer.resize(totalBytes);
        if (totalBytes > 0) memcpy(buffer.data(), flatIntBuffer.data(), totalBytes);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double bcastCommStartTime = MPI_Wtime();
    MPI_Bcast(&totalBytes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    if (totalBytes > 0) {
        if (rank != 0) buffer.resize(totalBytes);
        MPI_Bcast(buffer.data(), totalBytes, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    double bcastCommDuration = MPI_Wtime() - bcastCommStartTime;
    double maxBcastCommDuration; MPI_Reduce(&bcastCommDuration, &maxBcastCommDuration, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: bcastBlockNeighborMap_BcastComm_Max]: " << maxBcastCommDuration << "s\n";


    if (totalBytes == 0) {
        // std::cout << "Rank " << rank << ": broadcastBlockNeighborMap totalBytes is 0." << std::endl; // Moved
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Rank 0, Phase: bcastBlockNeighborMap_Total_Empty]: " << (MPI_Wtime() - totalPhaseStartTime) << "s\n";
        return (rank == 0) ? mapToSend : std::unordered_map<int, std::vector<int>>{};
    }

    std::unordered_map<int, std::vector<int>> receivedMap;
    if (rank != 0 && totalBytes > 0) {
        if (totalBytes < 0 || totalBytes % sizeof(int) != 0) MPI_Abort(MPI_COMM_WORLD,1);
        size_t numInts = static_cast<size_t>(totalBytes) / sizeof(int);
        std::vector<int> flatIntBuffer(numInts);
        if (totalBytes > 0) memcpy(flatIntBuffer.data(), buffer.data(), totalBytes);
        if (numInts == 0 && totalBytes > 0) MPI_Abort(MPI_COMM_WORLD,1);
        if (numInts > 0) {
            int numEntries = flatIntBuffer[0]; if (numEntries < 0) MPI_Abort(MPI_COMM_WORLD,1);
            size_t currentIndex = 1; receivedMap.reserve(numEntries);
            for (int i = 0; i < numEntries; ++i) {
                // Corrected Indentation
                if (currentIndex + 1 >= numInts) { 
                    std::cerr << "Rank " << rank << " Error: Buffer OOB reading key/count for entry " << i << ". Idx=" << currentIndex << ", Size=" << numInts << std::endl;
                    MPI_Abort(MPI_COMM_WORLD,1); 
                }
                int key = flatIntBuffer[currentIndex++];
                int numNeighbors = flatIntBuffer[currentIndex++];
                if (numNeighbors < 0) {
                    std::cerr << "Rank " << rank << " Error: Negative numNeighbors " << numNeighbors << " for key " << key << std::endl;
                    MPI_Abort(MPI_COMM_WORLD,1);
                }
                if (currentIndex + static_cast<size_t>(numNeighbors) > numInts) {
                    std::cerr << "Rank " << rank << " Error: Buffer overrun reading neighbors for key " << key << ". Need " << numNeighbors << ", Avail " << numInts - currentIndex << std::endl;
                    MPI_Abort(MPI_COMM_WORLD,1);
                }
                std::vector<int> neighbors;
                if (numNeighbors > 0) {
                    neighbors.assign(flatIntBuffer.begin() + currentIndex, flatIntBuffer.begin() + currentIndex + numNeighbors);
                    currentIndex += numNeighbors;
                }
                receivedMap[key] = std::move(neighbors);
            }
            if (currentIndex != numInts) {
                 std::cerr << "Rank " << rank << " Warning: Did not consume entire buffer in broadcastBlockNeighborMap. Index=" << currentIndex << ", Size=" << numInts << std::endl;
            }
        }
    }

    double totalDurationMap = MPI_Wtime() - totalPhaseStartTime;
    double maxTotalDurationMap; MPI_Reduce(&totalDurationMap, &maxTotalDurationMap, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) std::cout << std::fixed << std::setprecision(6) << "TIMING [Phase: bcastBlockNeighborMap_Total_Max]: " << maxTotalDurationMap << "s\n";
    // std::cout << "Rank " << rank << ": Received block neighbor map with " << ((rank == 0) ? mapToSend.size() : receivedMap.size()) << " entries.\n"; // Moved
    return (rank == 0) ? mapToSend : receivedMap;
}
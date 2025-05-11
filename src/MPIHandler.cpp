#include "../header/MPIHandler.h"
#include "../header/CSVParser.h" // Needed for mapToSIR
#include "../header/SIRCell.h"   // Needed for SIRCell type

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <list>
#include <numeric>        // For std::accumulate (maybe needed later)
#include <unordered_map>  // For broadcastBlockNeighborMap
#include <cstring>        // For memcpy
#include <stdexcept>      // For runtime_error (optional)
#include <set>            // For unique cell IDs in getDataForLocalBlocks

// Constructor: Initialize MPI
MPIHandler::MPIHandler(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

// Destructor: Finalize MPI
MPIHandler::~MPIHandler() {
    // Ensure MPI is finalized only if it was initialized successfully
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

// --- Method Implementations ---

// distributeData (original method, might be superseded by block distribution)
// Kept for completeness, but may not be used in the block-based workflow.
std::vector<SIRCell> MPIHandler::distributeData(const std::vector<std::vector<double>>& fullData) {
     std::vector<SIRCell> localGrid;
     int totalRows = 0;

     if (rank == 0) {
         totalRows = static_cast<int>(fullData.size());
         std::cout << "Rank 0: Distributing " << totalRows << " total rows (distributeData)." << std::endl;
     }
     // Broadcast total rows
     MPI_Bcast(&totalRows, 1, MPI_INT, 0, MPI_COMM_WORLD);

     if (totalRows == 0) {
         std::cout << "Rank " << rank << ": No data to distribute (distributeData)." << std::endl;
         return localGrid; // Return empty if no data
     }

     // Calculate distribution
     int rowsPerProc = totalRows / size;
     int extra = totalRows % size;
     int localRows = (rank < extra) ? (rowsPerProc + 1) : rowsPerProc;
     int startIndex = (rank < extra) ? (rank * (rowsPerProc + 1)) : (rank * rowsPerProc + extra);

     std::cout << "Rank " << rank << " assigned rows " << startIndex << " to "
               << (startIndex + localRows - 1) << " (" << localRows << " rows) (distributeData)." << std::endl;


     if (rank == 0) {
         // Rank 0 keeps its own data
         for (int i = startIndex; i < startIndex + localRows; ++i) {
              if (i >= 0 && static_cast<size_t>(i) < fullData.size()) {
                 try {
                    localGrid.push_back(CSVParser::mapToSIR(fullData[i]));
                 } catch (const std::exception& e) {
                     std::cerr << "Rank 0 Warning: mapToSIR failed for row " << i << " in distributeData: " << e.what() << std::endl;
                 }
              }
         }
         std::cout << "Rank 0 kept " << localGrid.size() << " rows (distributeData)." << std::endl;

         // Send data to other processes
         for (int proc = 1; proc < size; ++proc) {
             int procRows = (proc < extra) ? (rowsPerProc + 1) : rowsPerProc;
             int procStart = (proc < extra) ? (proc * (rowsPerProc + 1)) : (proc * rowsPerProc + extra);
             std::vector<double> sendBuffer;
             sendBuffer.reserve(procRows * 3); // S, I, R per row

             for (int i = procStart; i < procStart + procRows; ++i) {
                  if (i >= 0 && static_cast<size_t>(i) < fullData.size()) {
                     try {
                         SIRCell cell = CSVParser::mapToSIR(fullData[i]);
                         sendBuffer.push_back(cell.getS());
                         sendBuffer.push_back(cell.getI());
                         sendBuffer.push_back(cell.getR());
                     } catch (const std::exception& e) {
                          std::cerr << "Rank 0 Warning: mapToSIR failed for row " << i << " preparing send to rank " << proc << " in distributeData: " << e.what() << std::endl;
                     }
                  } else {
                      std::cerr << "Rank 0 Warning: Index " << i << " out of bounds for sending to rank " << proc << " (distributeData)." << std::endl;
                  }
             }

             // Send data size (number of doubles) first, then data
             int sendSize = sendBuffer.size();
             MPI_Send(&sendSize, 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
             if (sendSize > 0) {
                 MPI_Send(sendBuffer.data(), sendSize, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD);
                 std::cout << "Rank 0 sent " << sendSize / 3 << " rows (" << sendSize << " doubles) to rank " << proc << " (distributeData)." << std::endl;
             } else {
                  std::cout << "Rank 0 sent 0 rows to rank " << proc << " (distributeData)." << std::endl;
             }
         }
     } else {
         // Other processes receive data size, then data
         int recvSize;
         MPI_Recv(&recvSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

         if (recvSize > 0) {
             // Basic check for sanity of received size
             if (recvSize % 3 != 0) {
                  std::cerr << "Rank " << rank << " Error: Received size " << recvSize << " not divisible by 3 in distributeData. Aborting." << std::endl;
                  MPI_Abort(MPI_COMM_WORLD, 1);
             }
             std::vector<double> recvBuffer(recvSize);
             MPI_Recv(recvBuffer.data(), recvSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             std::cout << "Rank " << rank << " received " << recvSize / 3 << " rows (" << recvSize << " doubles) (distributeData)." << std::endl;

             // Unpack received data
             localGrid.reserve(recvSize / 3);
             for (int i = 0; i < recvSize / 3; ++i) {
                 try {
                     SIRCell cell(recvBuffer[3 * i], recvBuffer[3 * i + 1], recvBuffer[3 * i + 2]);
                     localGrid.push_back(cell);
                 } catch (const std::exception& e) { // Assuming SIRCell constructor could throw
                     std::cerr << "Rank " << rank << " Error constructing SIRCell from received data at index " << i << ": " << e.what() << std::endl;
                 }
             }
         } else {
              std::cout << "Rank " << rank << " received 0 rows (distributeData)." << std::endl;
         }
     }

     MPI_Barrier(MPI_COMM_WORLD); // Synchronize after distribution
     std::cout << "Rank " << rank << " finished distributeData with " << localGrid.size() << " local cells." << std::endl;
     return localGrid;
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
            std::cerr << "Rank " << rank << " Warning: Invalid step data size. Padding." << std::endl;
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

    // Normalize results on rank 0
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

// Modified writeResults: now receives correct counts/displacements
void MPIHandler::writeResults(const std::vector<double>& globalFlat,
                              const std::vector<int>& recvCounts,
                              const std::vector<int>& displacements) {
    if (rank == 0) {
        std::ofstream outfile("./data/simulation_results.csv");
        if (!outfile) {
            std::cerr << "Rank 0: Cannot open simulation_results.csv" << std::endl;
            return;
        }

        outfile << "Rank,Time,S_avg,I_avg,R_avg\n";
        int doublesPerStep = 4;
        for (int proc = 0; proc < size; ++proc) {
            int startIdx = displacements[proc];
            int numDoubles = recvCounts[proc];
            if (numDoubles % doublesPerStep != 0) {
                std::cerr << "Rank 0 Warning: Misaligned data size from rank " << proc << std::endl;
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
        std::cout << "Rank 0: Results written to simulation_results.csv" << std::endl;
    }
}

// distributeBlocks: Distributes ONLY block structure (IDs)
std::map<int, std::list<int>> MPIHandler::distributeBlocks(
    const std::map<int, std::list<int>>& allBlocks) {

     std::map<int, std::list<int>> localBlocks;
     int totalBlocks = 0;

     if (rank == 0) {
         totalBlocks = allBlocks.size();
         std::cout << "Rank 0: Distributing structure of " << totalBlocks << " blocks." << std::endl;
     }
     // Broadcast total number of blocks
     MPI_Bcast(&totalBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

     if (totalBlocks == 0) {
         std::cout << "Rank " << rank << ": No blocks to distribute." << std::endl;
         return localBlocks; // Nothing to distribute
     }

     // Determine which blocks this process gets based on contiguous assignment of block indices
     // This assumes an implicit ordering (e.g., iteration order on rank 0 or sorted IDs)
     int blocksPerProc = totalBlocks / size;
     int extraBlocks = totalBlocks % size;
     int numMyBlocks = (rank < extraBlocks) ? (blocksPerProc + 1) : blocksPerProc;
     int startBlockIndex = (rank < extraBlocks) ? (rank * (blocksPerProc + 1)) : (rank * blocksPerProc + extraBlocks);


     if (rank == 0) {
         // Assign own blocks based on index
         int currentBlockIndex = 0;
         for(const auto& [blockId, cellList] : allBlocks) { // Relies on map iteration order if not sorted
            if (currentBlockIndex >= startBlockIndex && currentBlockIndex < startBlockIndex + numMyBlocks) {
                localBlocks[blockId] = cellList; // Rank 0 keeps its own block data
            }
            currentBlockIndex++;
         }

         // Send block structures to other processes
         int blockCounter = 0;
         for (const auto& [blockId, cellList] : allBlocks) {
              // Determine target rank for this block based on its index in the iteration
              int targetRank = -1;
              int tempBlocksPerProc = totalBlocks / size;
              int tempExtraBlocks = totalBlocks % size;
              for(int r=0; r<size; ++r){
                  int rNumBlocks = (r < tempExtraBlocks) ? (tempBlocksPerProc + 1) : tempBlocksPerProc;
                  int rankStartBlockIndex = (r < tempExtraBlocks) ? (r * (tempBlocksPerProc + 1)) : (r * tempBlocksPerProc + tempExtraBlocks);
                  if(blockCounter >= rankStartBlockIndex && blockCounter < rankStartBlockIndex + rNumBlocks){
                      targetRank = r;
                      break;
                  }
              }

             if (targetRank > 0 && targetRank < size) { // Send only if target is another valid rank
                 // Serialize: blockId, numCells, cell1, cell2, ...
                 std::vector<int> blockData;
                 blockData.push_back(blockId);
                 blockData.push_back(static_cast<int>(cellList.size()));
                 blockData.insert(blockData.end(), cellList.begin(), cellList.end());

                 // Send size then data (using distinct tags)
                 int dataSize = blockData.size();
                 MPI_Send(&dataSize, 1, MPI_INT, targetRank, 0, MPI_COMM_WORLD); // Tag 0 for size
                 if (dataSize > 0) {
                     MPI_Send(blockData.data(), dataSize, MPI_INT, targetRank, 1, MPI_COMM_WORLD); // Tag 1 for data
                 }
             } else if (targetRank != 0) {
                  std::cerr << "Rank 0 Warning: Calculated invalid target rank " << targetRank << " for block index " << blockCounter << std::endl;
             }
             blockCounter++;
         }

         // Send termination signal (size = -1) to other ranks
         int terminateSignal = -1;
         for (int proc = 1; proc < size; ++proc) {
             MPI_Send(&terminateSignal, 1, MPI_INT, proc, 0, MPI_COMM_WORLD); // Use tag 0 (size tag)
         }

     } else { // Other ranks receive block structure
         while (true) {
             int dataSize;
             MPI_Status status;
             MPI_Recv(&dataSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); // Receive size with tag 0

             if (dataSize == -1) { // Termination signal
                 std::cout << "Rank " << rank << ": Received termination signal for blocks." << std::endl;
                 break;
             }

             if (dataSize > 0) {
                 std::vector<int> blockData(dataSize);
                 // Receive data with tag 1
                 MPI_Recv(blockData.data(), dataSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                 // Deserialize carefully
                 if (blockData.size() >= 2) { // Need at least blockId and numCells
                     int blockId = blockData[0];
                     int numCells = blockData[1];
                     // Check if received size matches expected size based on numCells
                     if (blockData.size() == static_cast<size_t>(2 + numCells)) {
                         std::list<int> cellList;
                         for (int i = 0; i < numCells; ++i) {
                             cellList.push_back(blockData[2 + i]);
                         }
                         localBlocks[blockId] = cellList; // Store received block
                     } else {
                          std::cerr << "Rank " << rank << " Error: Received block data size mismatch for block " << blockId
                                    << ". Expected " << 2+numCells << " ints, got " << blockData.size() << std::endl;
                     }
                 } else {
                      std::cerr << "Rank " << rank << " Error: Received block data too small (size " << blockData.size() << ")" << std::endl;
                 }
             } else if (dataSize == 0) {
                  std::cerr << "Rank " << rank << " Warning: Received block data size 0 from rank 0." << std::endl;
             }
             else { // dataSize < -1
                  std::cerr << "Rank " << rank << " Warning: Received invalid block data size " << dataSize << " from rank 0." << std::endl;
             }
         }
     }

     MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks have finished receiving blocks
     std::cout << "Rank " << rank << " finished distributeBlocks (structure only) with " << localBlocks.size() << " local blocks." << std::endl;
     return localBlocks;
}


// getDataForLocalBlocks: Fetches the necessary initial condition data rows
std::map<int, std::vector<double>> MPIHandler::getDataForLocalBlocks(
    const std::map<int, std::list<int>>& localBlocks,
    const std::vector<std::vector<double>>& fullData) { // fullData only valid/used on rank 0

    std::map<int, std::vector<double>> localCellData;

    // 1. Each process determines the UNIQUE cell IDs it needs data for
    std::set<int> neededCellIdsSet;
    for (const auto& [blockId, cellList] : localBlocks) {
        neededCellIdsSet.insert(cellList.begin(), cellList.end());
    }
    // Convert set to vector for sending (order matters for unpacking later)
    std::vector<int> neededCellIds(neededCellIdsSet.begin(), neededCellIdsSet.end());


    if (rank == 0) {
        // --- Rank 0: Manages requests and distributes data ---

        // Gather sizes of ID request vectors from all ranks
        int myRequestSize = neededCellIds.size();
        std::vector<int> requestSizes(size);
        MPI_Gather(&myRequestSize, 1, MPI_INT, requestSizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Calculate displacements for receiving all requests via GatherV
        std::vector<int> requestDispls(size);
        int totalRequestedIds = 0;
        requestDispls[0] = 0;
        totalRequestedIds = requestSizes[0];
        for (int i = 1; i < size; ++i) {
            requestDispls[i] = requestDispls[i - 1] + requestSizes[i - 1];
            totalRequestedIds += requestSizes[i];
        }
        std::cout << "Rank 0: Total unique cell IDs requested across all ranks: " << totalRequestedIds << std::endl;

        // Gather all requested cell IDs into one buffer on Rank 0
        std::vector<int> gatheredIdsBuffer(totalRequestedIds);
        // Rank 0 provides its own needed IDs for the GatherV operation
        MPI_Gatherv(neededCellIds.data(), myRequestSize, MPI_INT,
                    gatheredIdsBuffer.data(), requestSizes.data(), requestDispls.data(), MPI_INT,
                    0, MPI_COMM_WORLD);

        // Determine the number of doubles per cell from the actual data (if available)
        int doublesPerCell = 0;
        if (!fullData.empty() && !fullData[0].empty()) {
            doublesPerCell = fullData[0].size();
        } else if (totalRequestedIds > 0) {
             std::cerr << "Rank 0 Error: Cannot determine data size per cell (fullData is empty). Aborting." << std::endl;
             MPI_Abort(MPI_COMM_WORLD, 1);
        }
        std::cout << "Rank 0: Determined doubles per cell = " << doublesPerCell << std::endl;

        // Prepare data for Scatterv: Flatten all requested data rows
        std::vector<int> sendDataSizes(size);       // Size (num doubles) to send to each rank
        std::vector<int> sendDataDispls(size);      // Displacement in the flattened send buffer
        std::vector<double> flatSendDataBuffer;     // Buffer holding all data doubles to be scattered

        int currentSendDispl = 0;
        for (int targetRank = 0; targetRank < size; ++targetRank) {
            // Extract the requested IDs for this target rank
            int startIdx = requestDispls[targetRank];
            int numIds = requestSizes[targetRank];
            std::vector<double> rankDataBuffer; // Temp buffer for this rank's data doubles

            for (int k = 0; k < numIds; ++k) {
                int cellId = gatheredIdsBuffer[startIdx + k]; // Get the requested ID
                 if (cellId >= 0 && static_cast<size_t>(cellId) < fullData.size()) {
                     const auto& rowData = fullData[cellId];
                     // Ensure data consistency
                     if (static_cast<int>(rowData.size()) != doublesPerCell) {
                          std::cerr << "Rank 0 Error: Inconsistent row size for cell ID " << cellId << ". Expected "
                                    << doublesPerCell << ", got " << rowData.size() << ". Aborting." << std::endl;
                          MPI_Abort(MPI_COMM_WORLD, 1);
                     }
                     // Append the row's data to this rank's temporary buffer
                     rankDataBuffer.insert(rankDataBuffer.end(), rowData.begin(), rowData.end());

                     // If processing for rank 0 itself, also populate its localCellData map
                     if (targetRank == 0) {
                         localCellData[cellId] = rowData;
                     }
                 } else {
                     // Handle request for invalid/out-of-bounds cell ID
                     std::cerr << "Rank 0 Warning: Rank " << targetRank << " requested invalid cell ID " << cellId
                               << " (max index " << fullData.size() - 1 << "). Data will not be sent for this ID." << std::endl;
                     // We need to account for this missing data in sendDataSizes if we were packing ID+data
                     // Since we only pack data doubles now, the size mismatch will be handled by the receiver potentially.
                 }
            }
            // Store size and displacement for this rank's data segment
            sendDataSizes[targetRank] = rankDataBuffer.size(); // Total doubles for this rank
            sendDataDispls[targetRank] = currentSendDispl;
            // Append this rank's prepared data to the main send buffer
            flatSendDataBuffer.insert(flatSendDataBuffer.end(), rankDataBuffer.begin(), rankDataBuffer.end());
            currentSendDispl += rankDataBuffer.size(); // Update displacement for next rank
        } // End loop over target ranks


        // Broadcast the number of doubles per cell row so receivers know how to unpack
        MPI_Bcast(&doublesPerCell, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Broadcast the sizes array needed by Scatterv receivers
        MPI_Bcast(sendDataSizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

        // Scatter the prepared data segments to all ranks (including rank 0)
        std::vector<double> localRecvBuffer_Rank0(sendDataSizes[0]); // Rank 0 needs a buffer for its part
        MPI_Scatterv(flatSendDataBuffer.data(), sendDataSizes.data(), sendDataDispls.data(), MPI_DOUBLE,
                     localRecvBuffer_Rank0.data(), sendDataSizes[0], MPI_DOUBLE, // Rank 0 receives its segment
                     0, MPI_COMM_WORLD);
        // Note: Rank 0 already populated its localCellData map during the processing loop.

    } else { // --- Ranks > 0: Send request, receive data ---

        // Send size of needed IDs to rank 0
        int myRequestSize = neededCellIds.size();
        MPI_Gather(&myRequestSize, 1, MPI_INT, nullptr, 0, MPI_INT, 0, MPI_COMM_WORLD);

        // Send needed cell IDs to rank 0
        MPI_Gatherv(neededCellIds.data(), myRequestSize, MPI_INT,
                    nullptr, nullptr, nullptr, MPI_INT, // Non-root ranks don't receive the full buffer
                    0, MPI_COMM_WORLD);

        // Receive the number of doubles per cell row from rank 0
        int doublesPerCell = 0;
        MPI_Bcast(&doublesPerCell, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (doublesPerCell <= 0) {
             std::cerr << "Rank " << rank << " Error: Received invalid doublesPerCell (" << doublesPerCell << ") from rank 0. Aborting." << std::endl;
             MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Receive the array of data sizes (how much data each rank gets)
        std::vector<int> dataRecvSizes(size);
        MPI_Bcast(dataRecvSizes.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

        // Prepare local buffer to receive my segment of data
        int myRecvSize = dataRecvSizes[rank]; // Get my expected data size
        std::vector<double> localRecvBuffer(myRecvSize);

        // Receive my data segment using Scatterv
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, // Non-root ranks don't provide send buffer
                     localRecvBuffer.data(), myRecvSize, MPI_DOUBLE,
                     0, MPI_COMM_WORLD);

        // Unpack received data doubles into the localCellData map
        // Important: Assumes the order of doubles in localRecvBuffer corresponds
        // EXACTLY to the order of cell IDs in the 'neededCellIds' vector sent earlier.
        if (myRecvSize % doublesPerCell != 0 && myRecvSize != 0) {
             std::cerr << "Rank " << rank << " Error: Received data size " << myRecvSize
                       << " is not divisible by expected doubles per cell (" << doublesPerCell << "). Cannot unpack reliably." << std::endl;
             // Cannot proceed reliably if data size is wrong.
             MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
             size_t numCellsReceived = (myRecvSize == 0) ? 0 : static_cast<size_t>(myRecvSize / doublesPerCell);
             size_t bufferIdx = 0;

             // Check if the number of cells received matches the number requested (sanity check)
             if (neededCellIds.size() != numCellsReceived) {
                  std::cerr << "Rank " << rank << " Warning: Received data for " << numCellsReceived << " cells, but requested " << neededCellIds.size()
                            << ". This might happen if rank 0 skipped invalid requested IDs." << std::endl;
                   // Adjust loop limit based on received data to avoid buffer overruns
                  // This assumes Rank 0 packed data contiguously for the valid IDs it processed.
                  // The mapping back to cell IDs might be broken if IDs were skipped.
                  // A more robust solution involves rank 0 sending back which IDs it actually sent data for.
             }

             // Proceed assuming the order matches the request list, up to the number received
             for (size_t i = 0; i < numCellsReceived; ++i) {
                  // Check if we still have enough buffer space (should be guaranteed by size check above)
                  if (bufferIdx + doublesPerCell > static_cast<size_t>(myRecvSize)) {
                      std::cerr << "Rank " << rank << " Error: Buffer overrun during unpacking at cell index " << i << ". Aborting." << std::endl;
                      MPI_Abort(MPI_COMM_WORLD, 1);
                  }
                  // Get the cell ID corresponding to this data segment (assuming order)
                  int cellId = neededCellIds[i];
                  // Extract the row data for this cell
                  std::vector<double> rowData(localRecvBuffer.begin() + bufferIdx,
                                              localRecvBuffer.begin() + bufferIdx + doublesPerCell);
                  // Store in the local map
                  localCellData[cellId] = rowData;
                  bufferIdx += doublesPerCell; // Move buffer index to the next cell's data
             }
        }
    } // End else (rank > 0)

    MPI_Barrier(MPI_COMM_WORLD); // Ensure all ranks have finished receiving/processing
    std::cout << "Rank " << rank << " finished getDataForLocalBlocks. Found/received data for " << localCellData.size() << " cells." << std::endl;

    return localCellData;
}


// broadcastBlockNeighborMap: Broadcasts the block adjacency map
std::unordered_map<int, std::vector<int>> MPIHandler::broadcastBlockNeighborMap(
    const std::unordered_map<int, std::vector<int>>& mapToSend) {
    // ... Full implementation from previous step ...
    long long totalBytes = 0; std::vector<char> buffer;
    if (rank == 0) { /* Serialize map */
        std::vector<int> flatIntBuffer; flatIntBuffer.push_back(static_cast<int>(mapToSend.size()));
        for (const auto& [key, vec] : mapToSend) { flatIntBuffer.push_back(key); flatIntBuffer.push_back(static_cast<int>(vec.size())); flatIntBuffer.insert(flatIntBuffer.end(), vec.begin(), vec.end()); }
        totalBytes = static_cast<long long>(flatIntBuffer.size()) * sizeof(int); buffer.resize(totalBytes); if (totalBytes > 0) { memcpy(buffer.data(), flatIntBuffer.data(), totalBytes); }
        std::cout << "Rank 0: Broadcasting block neighbor map. Entries: " << mapToSend.size() << ", Total bytes: " << totalBytes << std::endl;
    }
    MPI_Bcast(&totalBytes, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    if (totalBytes == 0) { return (rank == 0) ? mapToSend : std::unordered_map<int, std::vector<int>>{}; }
    if (rank != 0) { buffer.resize(totalBytes); }
    if (totalBytes > 0) { MPI_Bcast(buffer.data(), totalBytes, MPI_BYTE, 0, MPI_COMM_WORLD); }
    std::unordered_map<int, std::vector<int>> receivedMap;
    if (rank != 0) { /* Deserialize map */
        if (totalBytes < 0 || totalBytes % sizeof(int) != 0) { std::cerr << "Rank " << rank << " Error: Invalid byte count " << totalBytes << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
        size_t numInts = static_cast<size_t>(totalBytes) / sizeof(int); std::vector<int> flatIntBuffer(numInts); if (totalBytes > 0) memcpy(flatIntBuffer.data(), buffer.data(), totalBytes);
        if (numInts == 0 && totalBytes > 0) {std::cerr<<"Rank "<<rank<<" Error: buffer copy failed?"<<std::endl; MPI_Abort(MPI_COMM_WORLD,1);}
        int numEntries = flatIntBuffer[0]; if (numEntries < 0) { std::cerr << "Rank " << rank << " Error: Received negative entries " << numEntries << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
        size_t currentIndex = 1; receivedMap.reserve(numEntries);
        for (int i = 0; i < numEntries; ++i) {
            if (currentIndex + 1 >= numInts) { std::cerr << "Rank " << rank << " Error: Buffer OOB reading key/count entry " << i << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            int key = flatIntBuffer[currentIndex++]; int numNeighbors = flatIntBuffer[currentIndex++]; if (numNeighbors < 0) { std::cerr << "Rank " << rank << " Error: Negative neighbors " << numNeighbors << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            if (currentIndex + numNeighbors > numInts) { std::cerr << "Rank " << rank << " Error: Buffer OOB reading neighbors key " << key << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); }
            std::vector<int> neighbors; if (numNeighbors > 0) { neighbors.assign(flatIntBuffer.begin() + currentIndex, flatIntBuffer.begin() + currentIndex + numNeighbors); currentIndex += numNeighbors; }
            receivedMap[key] = std::move(neighbors); // Use move
        }
        if (currentIndex != numInts) { std::cerr << "Rank " << rank << " Warning: Buffer not fully consumed. Idx=" << currentIndex << ", Size=" << numInts << std::endl; }
        std::cout << "Rank " << rank << " successfully received block neighbor map. Entries: " << receivedMap.size() << std::endl;
    }
    return (rank == 0) ? mapToSend : receivedMap;
}
#include "../../header/core/TimingUtils.h"  // Updated include path
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>   // For std::fixed, std::setprecision
#include <numeric>   // For std::accumulate
#include <algorithm> // For std::min_element, std::max_element
#include <fstream>   // For std::ofstream
#include <mpi.h>     // For MPI_Wtime etc.

namespace TimingUtils {

    // Static variable to hold the log file stream (only used by Rank 0)
    static std::ofstream timingLogFile;
    static bool logFileIsOpen = false;

    double startTimer() {
        // MPI_Barrier(MPI_COMM_WORLD); // Optional: synchronize before starting timer
        return MPI_Wtime();
    }

    double stopTimer(double startTime) {
        return MPI_Wtime() - startTime;
    }

    bool initLogFile(const std::string& filename) {
        // Only Rank 0 should open and manage the file
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            timingLogFile.open(filename);
            if (timingLogFile.is_open()) {
                logFileIsOpen = true;
                // Write header
                timingLogFile << "PhaseName,Statistic,Value,Units,NumRanks\n";
                std::cout << "TIMING_LOG: Initialized log file: " << filename << std::endl;
                return true;
            } else {
                std::cerr << "TIMING_LOG Error: Could not open log file: " << filename << std::endl;
                logFileIsOpen = false;
                return false;
            }
        }
        return true; // Other ranks report success as they don't manage the file
    }

    void logTimingSummary(const std::string& phaseName,
                          double minDuration,
                          double maxDuration,
                          double avgDuration,
                          int commSize) {
        if (logFileIsOpen) { // Assumes this is only called by Rank 0 after initLogFile
            timingLogFile << phaseName << ",Min," << minDuration << ",s," << commSize << "\n";
            timingLogFile << phaseName << ",Max," << maxDuration << ",s," << commSize << "\n";
            timingLogFile << phaseName << ",Avg," << avgDuration << ",s," << commSize << "\n";
            timingLogFile.flush(); // Ensure data is written
        }
    }
    
    void logRankSpecificTime(const std::string& phaseName,
                             double duration,
                             int rank) {
        if (logFileIsOpen) { // Assumes this is only called by Rank 0 after initLogFile
            timingLogFile << phaseName << ",Rank_" << rank << "_Time," << duration << ",s,1\n";
            timingLogFile.flush();
        }
    }


    void printAndLogTimingSummary(const std::string& phaseName,
                                  double localDuration,
                                  int currentRank,
                                  int commSize) {
        std::vector<double> allDurations;
        if (commSize <= 0) {
            if (currentRank == 0) {
                 std::cerr << "TIMING_SUMMARY Error [Phase: " << phaseName << "]: Invalid commSize (" << commSize << ")." << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            return;
        }

        if (currentRank == 0) {
            allDurations.resize(commSize);
        }

        MPI_Gather(&localDuration, 1, MPI_DOUBLE,
                   allDurations.data(), 1, MPI_DOUBLE,
                   0, MPI_COMM_WORLD);

        if (currentRank == 0) {
            if (!allDurations.empty()) {
                double minDuration = *std::min_element(allDurations.begin(), allDurations.end());
                double maxDuration = *std::max_element(allDurations.begin(), allDurations.end());
                double sumDuration = std::accumulate(allDurations.begin(), allDurations.end(), 0.0);
                double avgDuration = (commSize > 0) ? sumDuration / commSize : 0.0; // Added commSize > 0 check

                // Print to console
                std::cout << std::fixed << std::setprecision(6);
                std::cout << "TIMING_SUMMARY [Phase: " << phaseName << "]: "
                          << "Min: " << minDuration << "s, "
                          << "Max: " << maxDuration << "s, "
                          << "Avg: " << avgDuration << "s (across " << commSize << " ranks)" << std::endl;

                // Log to file
                logTimingSummary(phaseName, minDuration, maxDuration, avgDuration, commSize);
            } else {
                 std::cout << "TIMING_SUMMARY [Phase: " << phaseName << "]: No duration data gathered (allDurations empty)." << std::endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void printAndLogRank0PhaseTime(const std::string& phaseName,
                                   double phaseStartTime, // Start time recorded by Rank 0
                                   int currentRank) {
        if (currentRank == 0) {
            double duration = MPI_Wtime() - phaseStartTime; // End time is now
            // Print to console
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "TIMING [Rank 0, Phase: " << phaseName << "]: " << duration << "s" << std::endl;
            // Log to file
            logRankSpecificTime(phaseName, duration, 0);
        }
        // Calling function should handle barriers if needed for strict output ordering
        // or if other ranks depend on this timing being completed.
    }

    void printAndLogIndividualRankTime(const std::string& phaseName,
                                       double localDuration,
                                       int currentRank) {
        // Print to console
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "TIMING [Rank " << currentRank << ", Phase: " << phaseName << "]: " << localDuration << "s" << std::endl;

        // For logging individual rank times comprehensively, we'd typically gather them on rank 0.
        // For this simpler version, we can have Rank 0 log its own time if it calls this,
        // or modify it to send data to Rank 0 for logging.
        // Let's have it so if Rank 0 calls this for itself, it logs.
        if (currentRank == 0) {
            logRankSpecificTime(phaseName + "_Rank0", localDuration, 0);
        }
        // If you want ALL individual rank times logged, you'd need:
        // 1. Each rank sends its (phaseName, localDuration) to Rank 0.
        // 2. Rank 0 receives and calls logRankSpecificTime for each.
        // This is more complex than the current simple print.
    }


    void closeLogFile() {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0 && logFileIsOpen) {
            timingLogFile.close();
            logFileIsOpen = false;
            std::cout << "TIMING_LOG: Closed log file." << std::endl;
        }
    }

} // namespace TimingUtils
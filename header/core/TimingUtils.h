// header/TimingUtils.h
#ifndef TIMING_UTILS_H
#define TIMING_UTILS_H

#include <string>
#include <vector>
#include <fstream> // For std::ofstream
#include <mpi.h> // For MPI_Wtime, MPI_Gather, MPI_Barrier, MPI_COMM_WORLD, MPI_DOUBLE, MPI_MAX

namespace TimingUtils {

    // --- Core Timer Functions ---
    double startTimer();
    double stopTimer(double startTime);

    // --- Console Printing & File Logging Functions ---

    // Initializes the timing log file (called by Rank 0 at the beginning)
    // Returns true on success, false on failure to open file.
    bool initLogFile(const std::string& filename);

    // Logs a summary entry to the file (called by Rank 0)
    void logTimingSummary(const std::string& phaseName,
                          double minDuration,
                          double maxDuration,
                          double avgDuration,
                          int commSize);

    // Logs an individual or Rank 0 specific timing entry (called by Rank 0 if it's a Rank 0 specific time)
    void logRankSpecificTime(const std::string& phaseName,
                             double duration,
                             int rank);


    // Prints summary to console AND logs it (called by all, Rank 0 writes/prints)
    void printAndLogTimingSummary(const std::string& phaseName,
                                  double localDuration,
                                  int currentRank,
                                  int commSize);

    // Prints Rank 0 phase time to console AND logs it (called by all, Rank 0 writes/prints)
    void printAndLogRank0PhaseTime(const std::string& phaseName,
                                   double phaseStartTime, // Start time recorded by Rank 0
                                   int currentRank);

    // Prints individual rank time to console AND logs it (called by the specific rank, logged via Rank 0)
    // For individual rank times, logging might be better done by gathering to rank 0 first.
    // Let's simplify: This will print, and if currentRank is 0, it will also log its own time.
    // For more comprehensive individual logging, a gather mechanism would be needed.
    void printAndLogIndividualRankTime(const std::string& phaseName,
                                       double localDuration,
                                       int currentRank);


    // Closes the timing log file (called by Rank 0 at the end)
    void closeLogFile();


    // --- Deprecated or for console-only, if you want to keep them separate ---
    // void printTimingSummary(const std::string& phaseName, double localDuration, int currentRank, int commSize);
    // void printRank0PhaseTime(const std::string& phaseName, double phaseStartTime, int currentRank);
    // void printIndividualRankTime(const std::string& phaseName, double localDuration, int currentRank);


} // namespace TimingUtils

#endif // TIMING_UTILS_H
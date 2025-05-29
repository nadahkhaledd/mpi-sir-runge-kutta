#ifndef MPI_HANDLER_H
#define MPI_HANDLER_H

#include <vector>
#include <map>
#include <list>
#include <unordered_map>

class SIRCell;

class MPIHandler {
public:
    // Constructor / Destructor
    MPIHandler(int argc, char *argv[]);
    ~MPIHandler();

    // --- Getters ---
    int getRank() const;
    int getSize() const;

    // --- MPI Communication / Data Handling Methods ---
    std::vector<double> gatherResults(
        const std::vector<std::vector<double>>& localResults,
        std::vector<int>& recvCounts,
        std::vector<int>& displacements
    );

    void writeResults(
        const std::vector<double>& globalFlat,
        const std::vector<int>& recvCounts,
        const std::vector<int>& displacements
    );

    std::map<int, std::list<int>> distributeBlocks(
        const std::map<int, std::list<int>>& allBlocks
    );

    std::map<int, std::vector<double>> getDataForLocalBlocks(
        const std::map<int, std::list<int>>& localBlocks,
        const std::vector<std::vector<double>>& fullData
    );

    std::unordered_map<int, std::vector<int>> broadcastBlockNeighborMap(
        const std::unordered_map<int, std::vector<int>>& mapToSend
    );

private:
    int rank;
    int size;
    
};

#endif // MPI_HANDLER_H

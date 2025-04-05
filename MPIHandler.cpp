#include "MPIHandler.h"
#include "CSVParser.h"
#include <mpi.h>
#include <iostream>
#include <fstream>

MPIHandler::MPIHandler(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

MPIHandler::~MPIHandler() {
    MPI_Finalize();
}

int MPIHandler::getRank() const { 
    return rank; 
}

int MPIHandler::getSize() const { 
    return size; 
}

std::vector<SIRCell> MPIHandler::distributeData(const std::vector<std::vector<double>>& fullData) {
    std::vector<SIRCell> localGrid;
    
    int totalRows = 0;
    if (rank == 0) {
        totalRows = fullData.size();
    }
    MPI_Bcast(&totalRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int rowsPerProc = totalRows / size;
    int extra = totalRows % size;
    int startIndex, localRows;
    
    if (rank < extra) {
        localRows = rowsPerProc + 1;
        startIndex = rank * localRows;
    } else {
        localRows = rowsPerProc;
        startIndex = rank * localRows + extra;
    }
    
    if (rank == 0) {
        // Process 0 handles its own portion
        for (int i = startIndex; i < startIndex + localRows; i++) {
            localGrid.push_back(CSVParser::mapToSIR(fullData[i]));
        }
        
        // Send portions to other processes
        for (int proc = 1; proc < size; proc++) {
            int procRows = (proc < extra) ? rowsPerProc + 1 : rowsPerProc;
            int procStart = (proc < extra) ? proc * (rowsPerProc + 1) : proc * rowsPerProc + extra;
            std::vector<double> sendBuffer;
            
            for (int i = procStart; i < procStart + procRows; i++) {
                SIRCell cell = CSVParser::mapToSIR(fullData[i]);
                sendBuffer.push_back(cell.getS());
                sendBuffer.push_back(cell.getI());
                sendBuffer.push_back(cell.getR());
            }
            
            MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
        }
    } else {
        // Other processes receive their portion
        std::vector<double> recvBuffer(localRows * 3);
        MPI_Recv(recvBuffer.data(), localRows * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for (int i = 0; i < localRows; i++) {
            SIRCell cell(recvBuffer[3*i], recvBuffer[3*i+1], recvBuffer[3*i+2]);
            localGrid.push_back(cell);
        }
    }
    
    return localGrid;
}

std::vector<double> MPIHandler::gatherResults(const std::vector<std::vector<double>>& localResults) {
    int steps = localResults.size();
    std::vector<double> localFlat;
    
    for (auto &row : localResults) {
        localFlat.insert(localFlat.end(), row.begin(), row.end());
    }
    
    std::vector<double> globalFlat;
    if (rank == 0) {
        globalFlat.resize(steps * 4 * size);
    }
    
    MPI_Gather(localFlat.data(), steps * 4, MPI_DOUBLE,
               globalFlat.data(), steps * 4, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    
    return globalFlat;
}

void MPIHandler::writeResults(const std::vector<double>& globalFlat, int steps) {
    if (rank == 0) {
        std::ofstream outfile("simulation_results.csv");
        outfile << "Process,Time,S,I,R\n";
        
        for (int proc = 0; proc < size; proc++) {
            for (int i = 0; i < steps; i++) {
                int index = (proc * steps + i) * 4;
                outfile << proc << "," 
                        << globalFlat[index] << ","
                        << globalFlat[index+1] << ","
                        << globalFlat[index+2] << ","
                        << globalFlat[index+3] << "\n";
            }
        }
        
        outfile.close();
        std::cout << "Results written to simulation_results.csv" << std::endl;
    }
}
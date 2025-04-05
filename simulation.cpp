#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

// Model parameters
// More conservative model parameters
double beta = 0.3;      // Reduced from 0.5
double gammaRate = 0.1; 
double dt = 0.1;        // Reduced from 1.0 for numerical stability
int numSteps = 1000;

// SIR structure (normalized values, e.g., percentages)
struct SIR {
    double S, I, R;
};

// Helper function: trim whitespace from a string (if needed)
std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// Updated CSV parser that skips the header and parses the essential columns.
// Assumes CSV header is: Province_State,Country_Region,Last_Update,Lat,Long_,Confirmed,Deaths,Recovered,Active,Date,...
// and that columns of interest are in fixed positions.
std::vector<std::vector<double>> loadUSStateData(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file " << filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    std::string line;
    bool headerFound = false;
    int lineCount = 0;
    
    // Each row will be stored as a vector of double
    std::vector<std::vector<double>> data;
    
    // Process lines until we find valid data
    while (std::getline(infile, line)) {
        lineCount++;
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Check if this line appears to be the header
        if (line.find("Province_State") != std::string::npos) {
            std::cout << "Found header at line " << lineCount << ": " << line << std::endl;
            headerFound = true;
            continue;
        }
        
        // Skip any lines before the header
        if (!headerFound) {
            std::cout << "Skipping pre-header line: " << line << std::endl;
            continue;
        }
        
        // Process data line
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        // Split line on commas
        while (std::getline(ss, token, ',')) {
            tokens.push_back(trim(token));
        }
        
        // Verify we have enough columns
        if (tokens.size() < 9) {
            std::cerr << "Line " << lineCount << " has fewer than 9 columns: " << line << std::endl;
            continue;
        }
        
        try {
            // Convert values to doubles
            double lat = std::stod(tokens[3]);
            double lon = std::stod(tokens[4]);
            double confirmed = std::stod(tokens[5]);
            double deaths = std::stod(tokens[6]);
            double recovered = std::stod(tokens[7]);
            double active = std::stod(tokens[8]);
            
            // Store values
            data.push_back({lat, lon, confirmed, deaths, recovered, active});
        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid value at line " << lineCount << ": " << line << "\nError: " << e.what() << std::endl;
            continue;
        } catch (const std::out_of_range& e) {
            std::cerr << "Value out of range at line " << lineCount << ": " << line << "\nError: " << e.what() << std::endl;
            continue;
        }
    }
    
    std::cout << "Successfully parsed " << data.size() << " data rows from CSV." << std::endl;
    return data;
}
// Mapping dataset columns to SIR compartments.
// For simplicity, we assume that we normalize to 1 (i.e., percentages)
// Updated mapToSIR function
SIR mapToSIR(const std::vector<double>& rowData) {
    // rowData: [lat, lon, confirmed, deaths, recovered, active]  //first 3 rows ignored for numbered data
    
    // Calculate total population - if not available, estimate based on cases
    double totalPopulation = std::max(1000.0, rowData[2] + 1000.0); // Confirmed cases plus buffer
    
    // Active cases (I)
    double I = rowData[5] / totalPopulation; 
    
    // Recovered + deaths (R)
    double R = (rowData[3] + rowData[4]) / totalPopulation;
    
    // Susceptible (S) - ensure it's not negative
    double S = std::max(0.0, 1.0 - I - R);
    
    // Normalize to ensure S + I + R = 1
    double sum = S + I + R;
    if (sum > 0) {
        S /= sum;
        I /= sum;
        R /= sum;
    } else {
        // Fallback if all zeros
        S = 0.99;
        I = 0.01;
        R = 0.0;
    }
    
    // Ensure values are valid and within reasonable ranges
    S = std::max(0.0, std::min(1.0, S));
    I = std::max(0.0, std::min(1.0, I));
    R = std::max(0.0, std::min(1.0, R));
    
    return {S, I, R};
}
// RK4 step function for SIR cell dynamics
// Updated RK4 step function with stability checks
SIR rk4Step(const SIR &current) {
    // Ensure input values are valid
    SIR validCurrent = {
        std::max(0.0, std::min(1.0, current.S)),
        std::max(0.0, std::min(1.0, current.I)),
        std::max(0.0, std::min(1.0, current.R))
    };
    
    auto fS = [&](const SIR &state) -> double {
        return -beta * state.S * state.I;
    };
    auto fI = [&](const SIR &state) -> double {
        return beta * state.S * state.I - gammaRate * state.I;
    };
    auto fR = [&](const SIR &state) -> double {
        return gammaRate * state.I;
    };
    
    SIR k1, k2, k3, k4, next, temp;
    
    k1 = {dt * fS(validCurrent), dt * fI(validCurrent), dt * fR(validCurrent)};
    temp = {validCurrent.S + 0.5 * k1.S, validCurrent.I + 0.5 * k1.I, validCurrent.R + 0.5 * k1.R};
    k2 = {dt * fS(temp), dt * fI(temp), dt * fR(temp)};
    temp = {validCurrent.S + 0.5 * k2.S, validCurrent.I + 0.5 * k2.I, validCurrent.R + 0.5 * k2.R};
    k3 = {dt * fS(temp), dt * fI(temp), dt * fR(temp)};
    temp = {validCurrent.S + k3.S, validCurrent.I + k3.I, validCurrent.R + k3.R};
    k4 = {dt * fS(temp), dt * fI(temp), dt * fR(temp)};
    
    next = {validCurrent.S + (k1.S + 2*k2.S + 2*k3.S + k4.S) / 6.0,
            validCurrent.I + (k1.I + 2*k2.I + 2*k3.I + k4.I) / 6.0,
            validCurrent.R + (k1.R + 2*k2.R + 2*k3.R + k4.R) / 6.0};
    
    // Normalize to ensure S + I + R = 1
    double sum = next.S + next.I + next.R;
    if (sum > 0) {
        next.S /= sum;
        next.I /= sum;
        next.R /= sum;
    } else {
        // Fallback if all zeros
        next.S = 0.99;
        next.I = 0.01;
        next.R = 0.0;
    }
    
    // Final bounds check
    next.S = std::max(0.0, std::min(1.0, next.S));
    next.I = std::max(0.0, std::min(1.0, next.I));
    next.R = std::max(0.0, std::min(1.0, next.R));
    
    return next;
}

void updateGrid(std::vector<SIR> &grid) {
    std::vector<SIR> newGrid = grid;
    for (size_t i = 0; i < grid.size(); ++i) {
        newGrid[i] = rk4Step(grid[i]);
    }
    grid = newGrid;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Process 0 reads the entire dataset.
    std::vector<std::vector<double>> fullData;
    if (rank == 0) {
        fullData = loadUSStateData("initial_conditions.csv");
        std::cout << "Total rows in input dataset: " << fullData.size() << "\n";
    }
    
    // For this demonstration, assume we perform a simple grouping based on rank.
    // In practice, you would implement a clustering (or nearest-neighbor) algorithm using lat/long.
    // Here we simply partition the dataset rows among processes.
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
    
    // Each process creates its local vector of SIR values from its portion of the dataset.
    std::vector<SIR> localGrid;
    if (rank == 0) {
        // Process 0 scatters its own portion and then sends the rest to other processes.
        for (int i = startIndex; i < startIndex + localRows; i++) {
            localGrid.push_back(mapToSIR(fullData[i]));
        }
        // Now send portions to other processes
        for (int proc = 1; proc < size; proc++) {
            int procRows = (proc < extra) ? rowsPerProc + 1 : rowsPerProc;
            int procStart = (proc < extra) ? proc * (rowsPerProc + 1) : proc * rowsPerProc + extra;
            std::vector<double> sendBuffer;
            // Pack the SIR data as three doubles per row
            for (int i = procStart; i < procStart + procRows; i++) {
                SIR cell = mapToSIR(fullData[i]);
                sendBuffer.push_back(cell.S);
                sendBuffer.push_back(cell.I);
                sendBuffer.push_back(cell.R);
            }
            MPI_Send(sendBuffer.data(), sendBuffer.size(), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
        }
    } else {
        // Other processes receive their portion
        std::vector<double> recvBuffer(localRows * 3);
        MPI_Recv(recvBuffer.data(), localRows * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < localRows; i++) {
            SIR cell = {recvBuffer[3*i], recvBuffer[3*i+1], recvBuffer[3*i+2]};
            localGrid.push_back(cell);
        }
    }
    
    // Simulation: record evolution over time.
    // We will output time, S, I, R for each cell (or an average per block).
    // TODO: handle evolution of cases along each day since the initial data
    // Here, for simplicity, we update each cell independently.
    std::vector<std::vector<double>> localResults; // Each entry: [time, avg_S, avg_I, avg_R]
    for (int step = 0; step < numSteps; ++step) {
        // Update local grid
        updateGrid(localGrid);
        
        // Optionally: include spatial interactions among cells in the block
        
        // For demonstration, compute average S, I, R over the local block
        double sumS = 0, sumI = 0, sumR = 0;
        for (auto &cell : localGrid) {
            sumS += cell.S;
            sumI += cell.I;
            sumR += cell.R;
        }
        double avgS = sumS / localGrid.size();
        double avgI = sumI / localGrid.size();
        double avgR = sumR / localGrid.size();
        double timeVal = step * dt;
        localResults.push_back({timeVal, avgS, avgI, avgR});
        
        // Synchronize processes if needed
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // Gather results at process 0
    // (For simplicity, assume each process has the same number of time steps.)
    int steps = localResults.size();
    std::vector<double> localFlat;
    for (auto &row : localResults) {
        // Flatten each row [time, avg_S, avg_I, avg_R]
        localFlat.insert(localFlat.end(), row.begin(), row.end());
    }
    
    std::vector<double> globalFlat;
    if (rank == 0) {
        globalFlat.resize(steps * 4 * size);
    }
    MPI_Gather(localFlat.data(), steps * 4, MPI_DOUBLE,
               globalFlat.data(), steps * 4, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    
    // Process 0 writes the simulation evolution to a CSV file.
    if (rank == 0) {
        std::ofstream outfile("simulation_results.csv");
        outfile << "Process,Time,S,I,R\n";
        // Each process’s data is concatenated.
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
        std::cout << "Simulation evolution saved to simulation_results.csv\n";
    }
    
    MPI_Finalize();
    return 0;
}
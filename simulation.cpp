#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

// Model parameters (default values, can be read from CSV)
double beta = 0.5;   // transmission rate
double gammaRate = 0.1;  // recovery rate
double dt = 0.01;    // time step
int numSteps = 1000; // number of time steps

// Structure for storing SIR state per cell
struct SIR {
    double S;
    double I;
    double R;
};

// Function to parse a CSV file and load initial SIR data into a 2D grid
// Assumes CSV format: row, col, S, I, R (one row per cell)
std::vector<std::vector<SIR>> loadInitialData(const std::string& filename, int& nrows, int& ncols) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file " << filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string line;
    
    // Read the first line to get grid dimensions
    if (std::getline(infile, line)) {
        std::istringstream dimStream(line);
        if (!(dimStream >> nrows >> ncols)) {  // Ensure we correctly extract numbers
            std::cerr << "Invalid first line format in CSV (should be: nrows ncols)\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    } else {
        std::cerr << "Empty CSV file!\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Initialize grid with default values
    std::vector<std::vector<SIR>> grid(nrows, std::vector<SIR>(ncols, {1.0, 0.0, 0.0}));

    // Read each subsequent line
    while (std::getline(infile, line)) {
        std::istringstream ss(line);
        std::string token;
        int r, c;
        double S, I, R;

        // Extract and validate values
        if (!std::getline(ss, token, ',') || token.empty() || !(std::istringstream(token) >> r)) continue;
        if (!std::getline(ss, token, ',') || token.empty() || !(std::istringstream(token) >> c)) continue;
        if (!std::getline(ss, token, ',') || token.empty() || !(std::istringstream(token) >> S)) S = 1.0; // Default value
        if (!std::getline(ss, token, ',') || token.empty() || !(std::istringstream(token) >> I)) I = 0.0;
        if (!std::getline(ss, token, ',') || token.empty() || !(std::istringstream(token) >> R)) R = 0.0;

        // Ensure valid grid indices
        if (r >= 0 && r < nrows && c >= 0 && c < ncols) {
            grid[r][c] = {S, I, R};
        } else {
            std::cerr << "Warning: Skipping invalid row/column indices (" << r << "," << c << ")\n";
        }
    }

    return grid;
}


// RK4 step function for one cell (local SIR dynamics)
// Computes the new SIR state given current state, parameters and time step
SIR rk4Step(const SIR &current) {
    auto fS = [&](const SIR &state) -> double {
        return -beta * state.S * state.I;
    };
    auto fI = [&](const SIR &state) -> double {
        return beta * state.S * state.I - gammaRate * state.I;
    };
    auto fR = [&](const SIR &state) -> double {
        return gammaRate * state.I;
    };
    
    SIR k1, k2, k3, k4, next;
    
    // k1
    k1.S = dt * fS(current);
    k1.I = dt * fI(current);
    k1.R = dt * fR(current);
    
    // k2
    SIR temp;
    temp.S = current.S + 0.5 * k1.S;
    temp.I = current.I + 0.5 * k1.I;
    temp.R = current.R + 0.5 * k1.R;
    k2.S = dt * fS(temp);
    k2.I = dt * fI(temp);
    k2.R = dt * fR(temp);
    
    // k3
    temp.S = current.S + 0.5 * k2.S;
    temp.I = current.I + 0.5 * k2.I;
    temp.R = current.R + 0.5 * k2.R;
    k3.S = dt * fS(temp);
    k3.I = dt * fI(temp);
    k3.R = dt * fR(temp);
    
    // k4
    temp.S = current.S + k3.S;
    temp.I = current.I + k3.I;
    temp.R = current.R + k3.R;
    k4.S = dt * fS(temp);
    k4.I = dt * fI(temp);
    k4.R = dt * fR(temp);
    
    // Combine increments
    next.S = current.S + (k1.S + 2*k2.S + 2*k3.S + k4.S) / 6.0;
    next.I = current.I + (k1.I + 2*k2.I + 2*k3.I + k4.I) / 6.0;
    next.R = current.R + (k1.R + 2*k2.R + 2*k3.R + k4.R) / 6.0;
    
    return next;
}

// Update function for the entire grid (this can be parallelized with MPI)
// Here, we assume that each process handles a subgrid; boundary data exchange would be added.
void updateGrid(std::vector<std::vector<SIR>> &grid) {
    int nrows = grid.size();
    int ncols = grid[0].size();
    
    // Create a temporary grid to store new states
    std::vector<std::vector<SIR>> newGrid = grid;
    
    // For each cell, update using RK4 integration
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            // If needed, include spatial interactions with neighbors here.
            newGrid[i][j] = rk4Step(grid[i][j]);
        }
    }
    grid = newGrid;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // For simplicity, assume process 0 reads the CSV and broadcasts grid dimensions and data.
    int nrows = 0, ncols = 0;
    std::vector<std::vector<SIR>> grid;
    
    if (rank == 0) {
        // Provide the path to your CSV dataset (format: row,col,S,I,R)
        std::string csvFile = "initial_conditions.csv";
        grid = loadInitialData(csvFile, nrows, ncols);
        std::cout << "Loaded grid dimensions: " << nrows << " x " << ncols << "\n";
    }
    
    // Broadcast grid dimensions to all processes
    MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // For simplicity, assume we are dividing the grid equally among processes by rows.
    int rowsPerProc = nrows / size;
    int extra = nrows % size;
    int startRow, localRows;
    if (rank < extra) {
        localRows = rowsPerProc + 1;
        startRow = rank * localRows;
    } else {
        localRows = rowsPerProc;
        startRow = rank * localRows + extra;
    }
    
    // Each process extracts its subgrid
    std::vector<std::vector<SIR>> localGrid(localRows, std::vector<SIR>(ncols));
    if (rank == 0) {
        // Process 0 sends appropriate rows to other processes
        for (int proc = 1; proc < size; ++proc) {
            int procRows = (proc < extra) ? rowsPerProc + 1 : rowsPerProc;
            int procStart = (proc < extra) ? proc * (rowsPerProc + 1) : proc * rowsPerProc + extra;
            for (int i = 0; i < procRows; ++i) {
                // Flatten row data and send (each SIR has 3 doubles)
                MPI_Send(grid[procStart + i].data(), ncols * 3, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD);
            }
        }
        // Process 0 copies its own rows
        for (int i = 0; i < localRows; ++i) {
            localGrid[i] = grid[startRow + i];
        }
    } else {
        // Other processes receive their subgrid rows
        for (int i = 0; i < localRows; ++i) {
            MPI_Recv(localGrid[i].data(), ncols * 3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    
    // Main simulation loop: each process updates its local grid over numSteps time steps.
    for (int step = 0; step < numSteps; ++step) {
        // Optionally exchange ghost cell data with neighboring processes here
        
        updateGrid(localGrid);
        
        // Optionally, output intermediate results, synchronize, or gather data periodically.
        if (step % 100 == 0 && rank == 0) {
            std::cout << "Step " << step << " completed.\n";
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    // Gather local grids back to process 0 (not implemented fully here for brevity)
    // ...
    
    if (rank == 0) {
        // Process 0 can now output final simulation results (e.g., save to CSV)
        std::cout << "simulation is being saved\n";
        std::ofstream outfile("simulation_results.csv");
        outfile << "row,col,S,I,R\n";
        // For simplicity, assume grid was gathered into a global grid
        // Here we output from process 0's local grid as an example:
        for (int i = 0; i < localRows; ++i) {
            for (int j = 0; j < ncols; ++j) {
                outfile << (startRow + i) << "," << j << ","
                        << localGrid[i][j].S << "," 
                        << localGrid[i][j].I << "," 
                        << localGrid[i][j].R << "\n";
            }
        }
        std::cout << "simulation finished\n";
        outfile.close();
        std::cout << "Simulation results saved to simulation_results.csv\n";
    }
    
    MPI_Finalize();
    return 0;
}

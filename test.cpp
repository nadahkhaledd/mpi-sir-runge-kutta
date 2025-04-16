#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <limits>

// Define the State structure to hold geographic and epidemic data
struct State {
    std::string name;
    double lat, lon;
    double population;
    double susceptible;
    double infected;
    double recovered;
    int blockID;  // Which block this state belongs to
};

// SIR model structure
struct SIR {
    double S, I, R;
};

// Block structure for grouped states
struct Block {
    std::vector<int> stateIndices;  // Indices of states in this block
    double centerLat, centerLon;    // Center coordinates of the block
};

// Model parameters - can be adjusted
constexpr double BETA = 0.3;       // Infection rate
constexpr double GAMMA = 0.1;      // Recovery rate
constexpr double DT = 0.1;         // Time step
constexpr int NUM_STEPS = 100;     // Number of simulation steps
constexpr int TotalLenght = 10;     // Number of simulation steps
constexpr int NUM_BLOCKS = 10;     // Number of blocks to divide states into
constexpr double NEIGHBOR_INFLUENCE = 0.3;  // How much neighboring states influence infection
/*
// Function to calculate distance between two geographic points (in degrees)
double calcDistance(double lat1, double lon1, double lat2, double lon2) {
    return std::sqrt(std::pow(lat2 - lat1, 2) + std::pow(lon2 - lon1, 2));
}*/
//For better accuracy in calculating distances (especially considering the Earth's curvature), I suggested using the Haversine formula to calculate the distance between two points based on their latitude and longitude.
//The Haversine formula is more accurate for distances on a spherical Earth:

double calcDistance(double lat1, double lon1, double lat2, double lon2) {
    double R = 6371; // Earth's radius in km
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;

    lat1 = lat1 * M_PI / 180.0;
    lat2 = lat2 * M_PI / 180.0;

    double a = sin(dLat / 2) * sin(dLat / 2) +
               cos(lat1) * cos(lat2) * sin(dLon / 2) * sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    return R * c;
}
// Function to group states into blocks using a simple clustering approach
std::vector<Block> groupStatesIntoBlocks(std::vector<State>& states, int numBlocks) {
    std::vector<Block> blocks(numBlocks);
    
    // Initialize with random centers
    std::srand(std::time(nullptr));
    for (int i = 0; i < numBlocks; i++) {
        int randIdx = std::rand() % states.size();
        blocks[i].centerLat = states[randIdx].lat;
        blocks[i].centerLon = states[randIdx].lon;
    }

    bool changed = true;
    int maxIterations = 100;
    int iteration = 0;

    while (changed && iteration < maxIterations) {
        changed = false;
        iteration++;

        // Clear previous assignments
        for (auto& block : blocks) {
            block.stateIndices.clear();
        }

        // Assign each state to the nearest block
        for (int i = 0; i < states.size(); i++) {
            double minDist = std::numeric_limits<double>::max();
            int bestBlock = 0;

            for (int j = 0; j < numBlocks; j++) {
                double dist = calcDistance(states[i].lat, states[i].lon, 
                                           blocks[j].centerLat, blocks[j].centerLon);
                if (dist < minDist) {
                    minDist = dist;
                    bestBlock = j;
                }
            }

            blocks[bestBlock].stateIndices.push_back(i);
        }

        // Recalculate centers
        for (int i = 0; i < numBlocks; i++) {
            if (blocks[i].stateIndices.empty()) continue;

            double sumLat = 0, sumLon = 0;
            for (int stateIdx : blocks[i].stateIndices) {
                sumLat += states[stateIdx].lat;
                sumLon += states[stateIdx].lon;
            }

            double newCenterLat = sumLat / blocks[i].stateIndices.size();
            double newCenterLon = sumLon / blocks[i].stateIndices.size();

            if (std::abs(newCenterLat - blocks[i].centerLat) > 0.001 || 
                std::abs(newCenterLon - blocks[i].centerLon) > 0.001) {
                changed = true;
            }

            blocks[i].centerLat = newCenterLat;
            blocks[i].centerLon = newCenterLon;
        }
    }

    // ? Assign Block IDs to States
    for (int i = 0; i < numBlocks; i++) {
        for (int stateIdx : blocks[i].stateIndices) {
            states[stateIdx].blockID = i;
        }
    }

    // ? Debug: Print assigned block IDs
    std::cout << "Final State Assignments:\n";
    for (const auto& state : states) {
        std::cout << "State " << state.name << " assigned to Block " << state.blockID << "\n";
    }

    return blocks;
}

// Function to parse a CSV line into a State structure
bool parseCSVLine(const std::string& line, State& state) {
    std::istringstream ss(line);
    std::string field;
    std::vector<std::string> fields;
    
    // Parse CSV line into fields
    while (std::getline(ss, field, ',')) {
        fields.push_back(field);
    }
    
    // Check if we have enough fields
    if (fields.size() < 20) return false;
    
    try {
        state.name = fields[0]; // Province_State
        
        // Skip if not a US state
        if (fields[1] != "US") return false;
        
        // Parse lat/long
        if (!fields[3].empty() && !fields[4].empty()) {
            state.lat = std::stod(fields[3]);
            state.lon = std::stod(fields[4]);
        } else {
            return false; // Skip if no coordinates
        }
        
        // Parse cases
        if (!fields[5].empty() && !fields[6].empty() && !fields[8].empty()) {
            double confirmed = std::stod(fields[5]);
            double deaths = std::stod(fields[6]);
            double active = std::stod(fields[8]);
            double recovered = confirmed - deaths - active;
            
            // Ensure we have positive values
            if (confirmed < 0 || active < 0) return false;
            
            // Set SIR values (normalized to percentages)
            double totalPop = 100.0; // Using percentages
            state.infected = (active / totalPop);
            state.recovered = ((recovered + deaths) / totalPop);
            state.susceptible = 1.0 - state.infected - state.recovered;
            state.population = totalPop;
        } else {
            return false; // Skip if missing case data
        }
        
        return true;
    } catch (const std::exception& e) {
        return false; // Skip if any parsing errors
    }
}

// Function to load state data from CSV file
std::vector<State> loadStateData(const std::string& filename) {
    std::vector<State> states;
    std::ifstream infile(filename);
    
    if (!infile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    std::string line;
    // Skip header
    std::getline(infile, line);
    
    // Read data lines
    while (std::getline(infile, line)) {
        State state;
        if (parseCSVLine(line, state)) {
            states.push_back(state);
        }
    }
    
    return states;
}

// Function to calculate the influence of neighboring states on infection rates
double calculateNeighborInfluence(const std::vector<State>& states, int stateIdx,
                                 const std::vector<Block>& blocks) {
    double influence = 0.0;
    
    // Find which block this state belongs to
    int blockID = -1;
    for (int i = 0; i < blocks.size(); i++) {
        if (std::find(blocks[i].stateIndices.begin(), blocks[i].stateIndices.end(), stateIdx) 
            != blocks[i].stateIndices.end()) {
            blockID = i;
            break;
        }
    }
    
    if (blockID == -1) return influence;
    
    // Calculate influence from other states in the same block
    for (int neighborIdx : blocks[blockID].stateIndices) {
        if (neighborIdx == stateIdx) continue;
        
        double distance = calcDistance(states[stateIdx].lat, states[stateIdx].lon,
                                      states[neighborIdx].lat, states[neighborIdx].lon);
        
        // Add inverse-distance weighted influence
        double weight = 1.0 / (1.0 + distance);
        influence += weight * states[neighborIdx].infected;
    }
    
    return influence * NEIGHBOR_INFLUENCE;
}

// RK4 implementation for the SIR model with neighbor influence
SIR rk4Step(double S, double I, double R, double neighborInfluence,double deta,double gamma,double dt) {
    auto fS = [&](double s, double i) { return -deta * s * (i + neighborInfluence); };
    auto fI = [&](double s, double i) { return deta * s * (i + neighborInfluence) - gamma * i; };
    auto fR = [&](double i) { return gamma * i; };
    
    // First step
    double k1s = dt * fS(S, I);
    double k1i = dt * fI(S, I);
    double k1r = dt * fR(I);
    
    // Second step
    double k2s = dt * fS(S + 0.5 * k1s, I + 0.5 * k1i);
    double k2i = dt * fI(S + 0.5 * k1s, I + 0.5 * k1i);
    double k2r = dt * fR(I + 0.5 * k1i);
    
    // Third step
    double k3s = dt * fS(S + 0.5 * k2s, I + 0.5 * k2i);
    double k3i = dt * fI(S + 0.5 * k2s, I + 0.5 * k2i);
    double k3r = dt * fR(I + 0.5 * k2i);
    
    // Fourth step
    double k4s = dt * fS(S + k3s, I + k3i);
    double k4i = dt * fI(S + k3s, I + k3i);
    double k4r = dt * fR(I + k3i);
    
    // Result
    SIR result;
    result.S = S + (k1s + 2*k2s + 2*k3s + k4s) / 6.0;
    result.I = I + (k1i + 2*k2i + 2*k3i + k4i) / 6.0;
    result.R = R + (k1r + 2*k2r + 2*k3r + k4r) / 6.0;
    
    // Ensure bounds
    result.S = std::max(0.0, std::min(1.0, result.S));
    result.I = std::max(0.0, std::min(1.0, result.I));
    result.R = std::max(0.0, std::min(1.0, result.R));
    
    return result;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Master process loads data and distributes work
    std::vector<State> allStates;
    std::vector<State> OriginalStates;
    std::vector<Block> blocks;
    
    if (rank == 0) {
        std::cout << "Loading state data..." << std::endl;
        OriginalStates = loadStateData("covid_data.csv");
        allStates = OriginalStates;
        if (allStates.empty()) {
            std::cerr << "No valid state data found!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        std::cout << "Loaded " << allStates.size() << " states." << std::endl;
        
        // Group states into blocks
        blocks = groupStatesIntoBlocks(allStates, NUM_BLOCKS);
        
        // Print block assignments for verification
        for (int i = 0; i < blocks.size(); i++) {
            std::cout << "Block " << i << " contains " << blocks[i].stateIndices.size() 
                      << " states centered at (" << blocks[i].centerLat << ", " 
                      << blocks[i].centerLon << ")" << std::endl;
        }
    }
    
    // Broadcast the number of states
    int numStates = (rank == 0) ? allStates.size() : 0;
    MPI_Bcast(&numStates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // If not master, resize allStates
    if (rank != 0) {
        allStates.resize(numStates);
    }
    
    // Create an MPI datatype for State struct
    MPI_Datatype MPI_STATE;
    int blocklengths[7] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Aint displacements[7];
    MPI_Datatype types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    
    // We need to communicate only numerical data, not the name string
    MPI_Aint baseAddress;
    MPI_Get_address(&allStates[0].lat, &baseAddress);
    MPI_Get_address(&allStates[0].lat, &displacements[0]);
    MPI_Get_address(&allStates[0].lon, &displacements[1]);
    MPI_Get_address(&allStates[0].population, &displacements[2]);
    MPI_Get_address(&allStates[0].susceptible, &displacements[3]);
    MPI_Get_address(&allStates[0].infected, &displacements[4]);
    MPI_Get_address(&allStates[0].recovered, &displacements[5]);
    MPI_Get_address(&allStates[0].blockID, &displacements[6]);
    
    // Make relative to baseAddress
    for (int i = 0; i < 7; i++) {
        displacements[i] = MPI_Aint_diff(displacements[i], baseAddress);
    }
    
    MPI_Type_create_struct(7, blocklengths, displacements, types, &MPI_STATE);
    MPI_Type_commit(&MPI_STATE);
    
    // Broadcast state data (excluding string name)
    for (int i = 0; i < numStates; i++) {
        if (rank == 0) {
            // Master sends state data
            MPI_Bcast(&allStates[i].lat, 1, MPI_STATE, 0, MPI_COMM_WORLD);
        } else {
            // Other processes receive state data
            MPI_Bcast(&allStates[i].lat, 1, MPI_STATE, 0, MPI_COMM_WORLD);
        }
    }
    
    // Broadcast the number of blocks
    int numBlocks = NUM_BLOCKS;
    MPI_Bcast(&numBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Broadcast block info
    if (rank != 0) {
        blocks.resize(numBlocks);
    }
    
    for (int i = 0; i < numBlocks; i++) {
        // Broadcast block center coordinates
        MPI_Bcast(&blocks[i].centerLat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&blocks[i].centerLon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // Broadcast the number of states in this block
        int numStatesInBlock = (rank == 0) ? blocks[i].stateIndices.size() : 0;
        MPI_Bcast(&numStatesInBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // Resize and broadcast state indices
        if (rank != 0) {
            blocks[i].stateIndices.resize(numStatesInBlock);
        }
        MPI_Bcast(blocks[i].stateIndices.data(), numStatesInBlock, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    // Divide blocks among processes
    int blocksPerProcess = numBlocks / size;
    int extraBlocks = numBlocks % size;
    int startBlock = rank * blocksPerProcess + std::min(rank, extraBlocks);
    int endBlock = startBlock + blocksPerProcess + (rank < extraBlocks ? 1 : 0);
    
    std::cout << "Process " << rank << " handling blocks " << startBlock << " to " << (endBlock-1) << std::endl;
    
    // Simulation loop
    double beta_array[] = {0.4,0.2,0.1,0.05,0.01,0.005,0.001};
    double gamma_array[] = {0.4,0.2,0.1,0.05,0.01,0.005,0.001};
    double dt_array[] = {0.1,0.01,0.005,0.001};
    double Influence_array[] = {0.5,0.3,0.1,0.05,0.01};
for(int i=0;i<7;i++){
for(int j=0;j<7;j++){
for(int k=0;k<4;k++){
for(int m=0;m<5;m++){
double startTime = MPI_Wtime();
    for (int step = 0; step < TotalLenght/dt_array[k]; step++) {
    allStates = OriginalStates;
        // For each block assigned to this process
        for (int blockIdx = startBlock; blockIdx < endBlock; blockIdx++) {
            // Update each state in the block
            for (int stateIdx : blocks[blockIdx].stateIndices) {
                // Calculate influence from neighboring states
                double neighborInfluence = calculateNeighborInfluence(allStates, stateIdx, blocks);
                
                // Apply RK4 step
                SIR newValues = rk4Step(
                    allStates[stateIdx].susceptible,
                    allStates[stateIdx].infected,
                    allStates[stateIdx].recovered,
                    Influence_array[m], beta_array[i], gamma_array[j], dt_array[k]
                );
                
                // Update state values
                allStates[stateIdx].susceptible = newValues.S;
                allStates[stateIdx].infected = newValues.I;
                allStates[stateIdx].recovered = newValues.R;
            }
        }
        
        // All processes synchronize state values
        for (int i = 0; i < numStates; i++) {
            // Find which process is responsible for this state
            int ownerBlock = -1;
            int ownerProcess = -1;
            
            for (int b = 0; b < numBlocks; b++) {
                if (std::find(blocks[b].stateIndices.begin(), blocks[b].stateIndices.end(), i) 
                    != blocks[b].stateIndices.end()) {
                    ownerBlock = b;
                    break;
                }
            }
            
            if (ownerBlock != -1) {
                // Calculate which process owns this block
                int blocksPerProcess = numBlocks / size;
                int extraBlocks = numBlocks % size;
                for (int p = 0; p < size; p++) {
                    int pStartBlock = p * blocksPerProcess + std::min(p, extraBlocks);
                    int pEndBlock = pStartBlock + blocksPerProcess + (p < extraBlocks ? 1 : 0);
                    
                    if (ownerBlock >= pStartBlock && ownerBlock < pEndBlock) {
                        ownerProcess = p;
                        break;
                    }
                }
            }
            
            // Broadcast updated values from the owner process
            double values[3];
            if (rank == ownerProcess) {
                values[0] = allStates[i].susceptible;
                values[1] = allStates[i].infected;
                values[2] = allStates[i].recovered;
            }
            
            MPI_Bcast(values, 3, MPI_DOUBLE, ownerProcess, MPI_COMM_WORLD);
            
            if (rank != ownerProcess) {
                allStates[i].susceptible = values[0];
                allStates[i].infected = values[1];
                allStates[i].recovered = values[2];
            }
        }
        
        // Print progress every 10 steps
        /*if (rank == 0 && step % 10 == 0) {
            std::cout << "Step " << step << " completed." << std::endl;
        }*/
    }
    double endTime = MPI_Wtime();
    double totalI = 0;
  for (const auto& s : allStates) totalI += s.infected;
  double avgI = totalI / allStates.size();
    std::cout << "influence:"<<Influence_array[m]<< ",beta:"<<beta_array[i]<<",gamma:"<<gamma_array[j]<<",dt:"<<dt_array[k]
          << ",avgI:"<<avgI<<",time: " <<(endTime-startTime)<<std::endl;
}}}}
    // Output results
    if (rank == 0) {
        std::ofstream outfile("simulation_results.csv");
        outfile << "StateName,Lat,Long,Susceptible,Infected,Recovered,BlockID\n";
        
        for (int i = 0; i < allStates.size(); i++) {
            outfile << allStates[i].name << ","
                    << allStates[i].lat << ","
                    << allStates[i].lon << ","
                    << allStates[i].susceptible << ","
                    << allStates[i].infected << ","
                    << allStates[i].recovered << ","
                    << allStates[i].blockID << "\n";
        }
        
        outfile.close();
        std::cout << "Results written to simulation_results.csv" << std::endl;
    }
    
    // Free the MPI datatype
    MPI_Type_free(&MPI_STATE);
    
    MPI_Finalize();
    return 0;
}
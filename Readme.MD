# 🦠 SIR Simulation with MPI
This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI to simulate the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

## 📁 Project Structure

├── src/  
│   ├── SIRCell.cpp / .h         # Represents an individual cell in the grid (SIR logic per cell)  
│   ├── SIRModel.cpp / .h        # Manages the logic for the overall SIR simulation model  
│   ├── GridSimulation.cpp / .h  # Handles the 2D grid of cells and their interactions  
│   ├── MPIHandler.cpp / .h      # Manages MPI communication between processes  
│   ├── CSVParser.cpp / .h       # Parses input CSV files for initial conditions  
│   └── main.cpp                 # Entry point; sets up simulation and runs it  
├── scripts/  
│   └── sort_csv_by_states.py    # Script to preprocess and sort input CSV data by US states  
├── data/  
│   ├── initial_conditions.csv   # Input data for the simulation  
│   └── sorted_initial_conditions.csv # Preprocessed and sorted input data  
├── output/                      # Directory for simulation output  
├── Makefile                     # Build instructions  
└── README.md                    # This file  

## 🧠 File Descriptions

### src/SIRCell.cpp / src/SIRCell.h
Defines the behavior of a single cell in the simulation. Each cell holds:
- Susceptible, infected, and recovered counts
- Methods for updating state based on infection and recovery rates

### src/SIRModel.cpp / src/SIRModel.h
Encapsulates the simulation model logic:
- Defines infection and recovery rates
- Applies updates to all cells
- Handles time step progression

### src/GridSimulation.cpp / src/GridSimulation.h
Handles the 2D grid environment:
- Manages spatial relationships between cells
- Spreads infection between neighboring cells
- Iterates over time steps

### MPIHandler.cpp / MPIHandler.h
Abstracts away MPI communication:
- Divides the grid among processes
- Synchronizes boundary cells between neighboring processes
- Gathers results at the end of the simulation

### CSVParser.cpp / CSVParser.h
Handles input/output:
- Reads initial population and infection data from CSV files
- Optionally writes simulation output to CSV

### main.cpp
The main entry point:
- Initializes MPI
- Parses input arguments
- Sets up the simulation
- Runs the simulation loop
- Finalizes MPI

## 🛠️ Building the Project

### Steps to Compile:
1. Open a terminal and navigate to the project directory:
   ```bash
   cd /path/to/project
   ```
2. Clean any previous build artifacts:
   ```bash
   make clean
   ```
3. Compile the project:
   ```bash
   make
   ```

This will compile all `.cpp` files and produce the executable:

```bash
output/sir_simulation
```

## 🚀 Running the Simulation

Run the simulation using `mpirun` or `mpiexec`, specifying the number of processes. For example:

```bash
mpirun -np 4 ./sir_simulation [options]
```

### Notes:
- Ensure that MPI is installed on your system (e.g., OpenMPI or MPICH).
- Add any additional options as needed for your simulation.

## 🧹 Cleaning the Project
To clean up build artifacts (e.g., object files and the executable), run:
```bash
make clean
```
This will remove all files in the `output` directory and the compiled executable.

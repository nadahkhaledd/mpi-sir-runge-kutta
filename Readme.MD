# 🦠 SIR Simulation with MPI
This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI for simulating the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

## 📁 Project Structure

├── SIRCell.cpp / .h         # Represents an individual cell in the grid (SIR logic per cell)  
├── SIRModel.cpp / .h        # Manages the logic for the overall SIR simulation model  
├── GridSimulation.cpp / .h  # Handles the 2D grid of cells and their interactions  
├── MPIHandler.cpp / .h      # Manages MPI communication between processes  
├── CSVParser.cpp / .h       # Parses input CSV files for initial conditions  
├── main.cpp                 # Entry point; sets up simulation and runs it  
├── Makefile                 # Build instructions  
└── README.md                # This file  

## 🧠 File Descriptions

### SIRCell.cpp / SIRCell.h
Defines the behavior of a single cell in the simulation. Each cell holds:
- Susceptible, infected, recovered counts
- Methods for updating state based on infection and recovery rates

### SIRModel.cpp / SIRModel.h
Encapsulates the simulation model logic:
- Defines infection and recovery rates
- Applies updates to all cells
- Handles time step progression

### GridSimulation.cpp / GridSimulation.h
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
mpirun -np 4 ./output/sir_simulation ./data/<input_file.csv> [options]
```

### Example Command:
```bash
mpirun -np 4 ./output/sir_simulation ./data/initial_conditions.csv --infection-rate 0.1 --recovery-rate 0.05 --steps 100
```

### Arguments:
- `<input_file.csv>`: Input CSV file with initial conditions.
- `--infection-rate`: Infection rate (e.g., `0.1`).
- `--recovery-rate`: Recovery rate (e.g., `0.05`).
- `--steps`: Number of simulation steps (e.g., `100`).

### Notes:
- Ensure that MPI is installed on your system (e.g., OpenMPI or MPICH).
- Replace `<input_file.csv>` with the path to your input file.
- Add any additional options as needed for your simulation.

## 🧹 Cleaning the Project
To clean up build artifacts (e.g., object files and the executable), run:
```bash
make clean
```
This will remove all files in the `output` directory and the compiled executable.

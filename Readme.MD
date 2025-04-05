🦠 # SIR Simulation with MPI
This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI for simulating the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

📁 ## Project Structure

`
├── SIRCell.cpp / .h         # Represents an individual cell in the grid (SIR logic per cell)

├── SIRModel.cpp / .h        # Manages the logic for the overall SIR simulation model

├── GridSimulation.cpp / .h  # Handles the 2D grid of cells and their interactions

├── MPIHandler.cpp / .h      # Manages MPI communication between processes

├── CSVParser.cpp / .h       # Parses input CSV files for initial conditions

├── main.cpp                 # Entry point; sets up simulation and runs it

├── Makefile                 # Build instructions

└── README.md                # This file
`
## 🧠 File Descriptions
SIRCell.cpp / SIRCell.h
Defines the behavior of a single cell in the simulation. Each cell holds:

Susceptible, infected, recovered counts
Methods for updating state based on infection and recovery rates
SIRModel.cpp / SIRModel.h
Encapsulates the simulation model logic:

Defines infection and recovery rates
Applies updates to all cells
Handles time step progression
GridSimulation.cpp / GridSimulation.h
Handles the 2D grid environment:

Manages spatial relationships between cells
Spreads infection between neighboring cells
Iterates over time steps
MPIHandler.cpp / MPIHandler.h
Abstracts away MPI communication:

Divides the grid among processes
Synchronizes boundary cells between neighboring processes
Gathers results at the end of the simulation
CSVParser.cpp / CSVParser.h
Handles input/output:

Reads initial population and infection data from CSV files
Optionally writes simulation output to CSV
main.cpp
The main entry point:

Initializes MPI
Parses input arguments
Sets up the simulation
Runs the simulation loop
Finalizes MPI

 ## 🛠️ Building the Project
Make sure you have MPI installed (e.g., OpenMPI or MPICH).


make
This will compile all .cpp files and produce the executable:


sir_simulation

## 🚀 Running the Simulation
Run with mpirun or mpiexec, specifying the number of processes:


mpirun -np 4 ./sir_simulation <input_file.csv> [options]
Arguments may include:

Input CSV file with initial conditions
Simulation parameters such as infection rate, recovery rate, number of steps, etc.

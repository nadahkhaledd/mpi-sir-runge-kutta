# 🦠 SIR Simulation with MPI
This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI for simulating the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

## 📁Project Structure

├── SIRCell.cpp / .h         # Represents an individual cell in the grid (SIR logic per cell)

├── SIRModel.cpp / .h        # Manages the logic for the overall SIR simulation model

├── GridSimulation.cpp / .h  # Handles the 2D grid of cells and their interactions

├── MPIHandler.cpp / .h      # Manages MPI communication between processes

├── CSVParser.cpp / .h       # Parses input CSV files for initial conditions

├── main.cpp                 # Entry point; sets up simulation and runs it

├── Makefile                 # Build instructions

├── test.cpp                 # Use Grid search method to findout the optimal hyperparameter(gamma beta..)

├── testResult.csv           # Show the time cost and accuracy by different set of parameter

├── covid_data.csv           # dataset for test.cpp to find out the hyperparameter

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


The testResult.csv record the cost(running time) and accuracy(eval based on averageInfectious)

The test result report shows the accuracy and code running efficiency of different parameters under different value combinations.

In general, according to the test results and data analysis, the change of dt will not cause huge cumulative errors in the accuracy of the results, but will significantly increase the computing cost; gamma and beta will not have too much impact on the computing cost, but when these two parameters are large, it will cause numerical instability and lead to numerical errors.


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

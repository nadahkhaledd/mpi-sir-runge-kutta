# 🦠 SIR Simulation with MPI
This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI to simulate the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

---

## 📁 Project Structure

├── src/  
│   ├── SIRCell.cpp / .h         # Represents an individual cell in the grid (SIR logic per cell)  
│   ├── SIRModel.cpp / .h        # Manages the logic for the overall SIR simulation model  
│   ├── GridSimulation.cpp / .h  # Handles the 2D grid of cells and their interactions  
│   ├── MPIHandler.cpp / .h      # Manages MPI communication between processes  
│   ├── CSVParser.cpp / .h       # Parses input CSV files for initial conditions  
│   └── main.cpp                 # Entry point; sets up simulation and runs it  
├── scripts/  
│   ├── sort_csv_by_states.py    # Script to preprocess and sort input CSV data by US states  
│   └── PlottingSIRModelResults.py # Script to plot the simulation results  
├── data/  
│   ├── initial_conditions.csv   # Input data for the simulation  
│   └── sorted_initial_conditions.csv # Preprocessed and sorted input data  
├── output/                      # Directory for simulation output  
├── Makefile                     # Build instructions  
└── README.md                    # This file

---

## 🧠 Key Tasks and How They Are Done

### 1️⃣ Dividing the Data
The simulation divides the grid of cells into blocks, which are distributed among MPI processes. The division is based on:
- **Number of processes**: The grid is divided into a number of blocks proportional to the number of MPI processes.
- **Optimal block configuration**: The `GridSimulation::divideIntoOptimalBlocks` method calculates the best block distribution by evaluating configurations based on:
  - Grid shape (how "square" the grid is).
  - Block shape (how "square" each block is).
  - Balance between the number of blocks and the size of each block.
  - Penalties for extreme configurations (e.g., too many blocks or too few cells per block).

Each block contains a list of cell IDs, and the cells are assigned to blocks in a way that ensures an even distribution of workload.

---

### 2️⃣ MPI Communication: Distributing, Broadcasting, and Gathering

#### **Distributing Blocks**
- The `MPIHandler::distributeBlocks` method divides the blocks among MPI processes.
- Rank 0 (the root process) determines the block distribution and sends the relevant blocks to each process using `MPI_Send`.
- Non-root processes receive their assigned blocks using `MPI_Recv`.

#### **Broadcasting Neighbor Maps**
- The `MPIHandler::broadcastBlockNeighborMap` method broadcasts the block adjacency map (neighbor relationships between blocks) from Rank 0 to all other ranks using `MPI_Bcast`.
- This ensures that all processes have the same view of the block relationships.

#### **Gathering Results**
- The `MPIHandler::gatherResults` method collects simulation results from all processes.
- Each process sends its local results (e.g., time step data for its cells) to Rank 0 using `MPI_Gatherv`.
- Rank 0 aggregates the results and normalizes the data (e.g., ensuring that the sum of S, I, and R values equals 1 for each cell).

---

### 3️⃣ Rank Correlation with Code Implementation

- **Rank 0**:
  - Responsible for loading the initial data from the CSV file.
  - Divides the grid into blocks and determines the optimal block distribution.
  - Broadcasts the block neighbor map to all other ranks.
  - Gathers and writes the final results to a CSV file.

- **Other Ranks**:
  - Receive their assigned blocks and the corresponding cell data.
  - Perform the simulation for their local cells.
  - Exchange boundary data with neighboring ranks to ensure consistency at block boundaries.

The rank of each process is accessed using the `MPIHandler::getRank` method, and the total number of processes is accessed using `MPIHandler::getSize`.

---

### 4️⃣ How Everything Is Put Together

1. **Initialization**:
   - MPI is initialized in `main.cpp` using the `MPIHandler` class.
   - Rank 0 loads the initial condition data from a CSV file.

2. **Data Division**:
   - The grid is divided into blocks using `GridSimulation::divideIntoOptimalBlocks`.
   - Blocks are distributed among MPI processes using `MPIHandler::distributeBlocks`.

3. **Simulation Setup**:
   - Each process receives its assigned blocks and the corresponding cell data.
   - Neighbor maps are built to track relationships between blocks and cells.
   - The grid is initialized with the local data for each process.

4. **Simulation Loop**:
   - Each process runs the simulation for its local cells using `GridSimulation::runSimulation`.
   - Boundary data is exchanged between neighboring processes to ensure consistency.

5. **Result Gathering**:
   - Each process sends its local results to Rank 0 using `MPIHandler::gatherResults`.
   - Rank 0 aggregates the results and writes them to a CSV file.

6. **Finalization**:
   - MPI is finalized, and the simulation ends.

---

### 5️⃣ Scoring Logic for Dividing Data

The `GridSimulation::divideIntoOptimalBlocks` method evaluates multiple configurations for dividing the grid into blocks and selects the best one based on a scoring system. Here's how the scoring works:

#### **1. Configurations Evaluated**
- The method considers all divisors of the total number of cells as potential numbers of blocks.
- For each divisor (number of blocks), it calculates the number of cells per block.

#### **2. Factors Considered in Scoring**
Each configuration is scored based on the following factors:

1. **Grid Shape**:
   - The overall grid of blocks is evaluated for how "square" it is.
   - A grid with a ratio close to 1.0 (e.g., 5x5) is preferred over a long, narrow grid (e.g., 1x25).

2. **Block Shape**:
   - Each block is evaluated for how "square" it is.
   - A block with a ratio close to 1.0 (e.g., 5x5) is preferred over a long, narrow block (e.g., 1x25).

3. **Cell Block Shape**:
   - The shape of the cells within each block is evaluated for how "square" it is.
   - A block with cells arranged in a square-like shape is preferred.

4. **Balance Factor**:
   - The balance between the number of blocks and the size of each block is evaluated.
   - Configurations with a reasonable balance (e.g., 10 blocks with 5 cells each) are preferred over extremes (e.g., 1 block with 50 cells or 50 blocks with 1 cell each).

5. **Penalties for Extreme Configurations**:
   - Configurations with too few blocks or too few cells per block are penalized.
   - For example, configurations with fewer blocks than the number of processes or fewer than 5 cells per block are penalized.

#### **3. Scoring Formula**
The score for each configuration is calculated as:
```
score = (gridRatio * 0.25) + (blockRatio * 0.25) + (cellRatio * 0.2) + (balanceFactor * 0.2) - penalty
```
Where:
- `gridRatio`: How square the overall grid of blocks is.
- `blockRatio`: How square each block is.
- `cellRatio`: How square the cells within each block are.
- `balanceFactor`: The balance between the number of blocks and the size of each block.
- `penalty`: A penalty applied for extreme configurations.

#### **4. Example Output**
For 50 cells and 4 processes, the method might produce the following output:
```
Finding optimal block distribution for 50 cells...
  Option: 1 blocks with 50 cells each. Grid: 5x10 Score: 0.45
  Option: 2 blocks with 25 cells each. Grid: 5x10 Score: 0.62
  Option: 5 blocks with 10 cells each. Grid: 2x25 Score: 0.68
  Option: 10 blocks with 5 cells each. Grid: 2x25 Score: 0.72
  Option: 25 blocks with 2 cells each. Grid: 5x10 Score: 0.55
  Option: 50 blocks with 1 cells each. Grid: 5x10 Score: 0.48
Optimal distribution: 10 blocks with 5 cells each.
```

#### **5. Selected Configuration**
The configuration with the highest score is selected as the optimal distribution. In this example, the configuration with 10 blocks and 5 cells per block is chosen.

---

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
sir_simulation
```

---

## 🚀 Running the Simulation

Run the simulation using `mpirun` or `mpiexec`, specifying the number of processes. For example:

```bash
mpirun -np 4 ./sir_simulation
```

### Notes:
- Ensure that MPI is installed on your system (e.g., OpenMPI or MPICH).
- Add any additional options as needed for your simulation.

---

## 📊 Plotting the Results

After running the simulation, you can visualize the results using the provided Python script. Follow these steps:

1. Navigate to the `scripts` directory:
   ```bash
   cd scripts
   ```

2. Run the plotting script:
   ```bash
   python PlottingSIRModelResults.py
   ```

### Requirements:
Ensure you have Python installed along with the required libraries (e.g., `matplotlib`, `numpy`). You can install the dependencies using:
```bash
pip install -r requirements.txt
```

---

## 🧹 Cleaning the Project
To clean up build artifacts (e.g., object files and the executable), run:
```bash
make clean
```
This will remove all files in the `output` directory and the compiled executable.



--

## Acknowledgements

We are deeply grateful to **Professor [Serena Curzel](https://www.linkedin.com/in/serena-curzel/)** and **Professor [Marco Gribaudo](https://www.linkedin.com/in/marco-gribaudo-41a327b3/)** for proposing the project idea and offering their continuous support and guidance throughout its development.

This hands-on challenge was part of the final exam for the course Advanced Methods for Scientific Computing in High-Performance Computing (2024/25) at Politecnico di Milano.

We also thank **Professor [Luca Formaggia](https://formaggia.faculty.polimi.it/)** for his exceptional teaching and mentorship throughout the course.


---

## Contributors

- [Salvatore Mariano Librici](https://www.linkedin.com/in/salvatore-mariano-librici/)
- [Nadah Khaled](https://www.linkedin.com/in/nadahkhaledd10/)
- [Milica Sanjevic](https://www.linkedin.com/in/milica-sanjevic-321392327/)
- [Hirdesh Kumar](https://www.linkedin.com/in/hirdeshkumar2407/)
- [Yibo Li](https://www.linkedin.com/in/yibo-li-0b0b792a2/)

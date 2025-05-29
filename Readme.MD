# 🦠 SIR Simulation with MPI

This project implements a parallel SIR (Susceptible-Infected-Recovered) model using MPI to simulate the spread of infectious diseases across a grid of cells. Each cell represents a population unit, and the simulation evolves over discrete time steps.

## Prerequisites & Installation

### Required Libraries
1. MPI Implementation
   ```bash
   # Ubuntu/Debian
   sudo apt-get update
   sudo apt-get install openmpi-bin libopenmpi-dev

   # CentOS/RHEL
   sudo yum install openmpi openmpi-devel
   ```

2. C++ Compiler (C++17 or later)
   ```bash
   sudo apt-get install g++
   ```

3. Python 3.x & Libraries
   ```bash
   sudo apt-get install python3 python3-pip
   pip3 install pandas matplotlib numpy
   ```

## 📁 Project Structure
```
├── src/  
│   ├── core/
│   │   ├── SIRCell.cpp/.h        # Individual cell in the grid
│   │   ├── SIRModel.cpp/.h       # Core SIR simulation model
│   │   ├── GridSimulation.cpp/.h # 2D grid and cell interactions
│   │   ├── MPIHandler.cpp/.h     # MPI communication management
│   │   ├── CSVParser.cpp/.h      # Input data parsing
│   │   └── SimulationManager.cpp/.h # High-level simulation control
│   └── tests/                    # Test implementations
├── scripts/  
│   ├── clean_sort_dataset.py     # Data preprocessing
│   └── PlottingSIRModelResults.py # Result visualization
├── data/  
│   ├── test_datasets/           # Test input files
│   │   └── README.md            # Test data documentation
│   ├── test_results/           # Test output files
│   └── output/                 # Main simulation output
└── Readme.MD                   # This file
```

## Building & Running

### Main Simulation
```bash
# Build
make clean
make simulation

# Run (e.g., with 4 processes)
mpirun -np 4 ./sir_simulation
```

### Test Suite
```bash
# Build tests
make clean
make test

# Run tests (e.g., with 4 processes)
mpirun -np 4 ./sir_test_suite
```

### Data Preprocessing
```bash
# Clean and sort a dataset
python scripts/clean_sort_dataset.py ./data/test_datasets/input.csv
```

## Testing Framework

### Available Tests
1. Temporal Tests
   - (Jan 2021)
   - (Feb 2021)
   - Tests simulation behavior across different time periods

2. Parameter Sensitivity Tests
   - Low beta (β=0.1)
   - High beta (β=0.5)
   - Tests model response to different transmission rates

### Running Tests
```bash
# Build test suite
make clean
make test

# Run all tests with 4 processes
mpirun -np 4 ./sir_test_suite

# Run with different number of processes
mpirun -np <num_processes> ./sir_test_suite
```

### Test Results
- Output Location: `./data/test_results/`
- Files Generated:
  - `<test_name>_p<num_processes>_results.csv`: Simulation results (e.g., `jan_2021_p4_results.csv` for a test run using 4 MPI processes)
  - `timing_log.csv`: Performance measurements for each process

### Adding New Tests
1. Add test dataset to `./data/test_datasets/`
2. Process dataset:
   ```bash
   python scripts/clean_sort_dataset.py ./data/test_datasets/new_dataset.csv
   ```
3. Add test configuration in `src/tests/TestSuite.cpp`:
   ```cpp
   configs.emplace_back(
       dataDir + "sorted_new_dataset.csv",
       resultDir + "new_test_name",
       beta, gamma, dt, steps, minProcs, maxProcs
   );
   ```

### Test Dataset Requirements
- Must follow format specified in data/test_datasets/README.md
- Include population and infection data for all 50 US states
- Headers must match expected format
- First line must contain dimensions: "<num_rows> <num_columns>"

## Input/Output Formats

### Input Dataset Requirements
- First line: "<num_rows> <num_columns>"
- Second line: Column headers
- Required columns:
  - Province_State
  - Population
  - Last_Update
  - Lat, Long_
  - Confirmed, Deaths, Recovered, Active

### Output Files
- Main simulation: `./data/output/simulation_results.csv`
- Test results: `./data/test_results/<test_name>_p<num_processes>_results.csv`
- Timing logs: 
  - Main: `./data/output/timing_log.csv`
  - Tests: `./data/test_results/timing_log.csv`

## Implementation Details

### 🧠 Key Tasks and How They Are Done

#### 1️⃣ Dividing the Data
The simulation divides the grid of cells into blocks, which are distributed among MPI processes. The division is based on:
- **Number of processes**: The grid is divided into a number of blocks proportional to the number of MPI processes.
- **Optimal block configuration**: The `GridSimulation::divideIntoOptimalBlocks` method calculates the best block distribution by evaluating configurations based on:
  - Grid shape (how "square" the grid is).
  - Block shape (how "square" each block is).
  - Balance between the number of blocks and the size of each block.
  - Penalties for extreme configurations (e.g., too many blocks or too few cells per block).

Each block contains a list of cell IDs, and the cells are assigned to blocks in a way that ensures an even distribution of workload.

---

#### 2️⃣ MPI Communication: Distributing, Broadcasting, and Gathering

##### **Distributing Blocks**
- The `MPIHandler::distributeBlocks` method divides the blocks among MPI processes.
- Rank 0 (the root process) determines the block distribution and sends the relevant blocks to each process using `MPI_Send`.
- Non-root processes receive their assigned blocks using `MPI_Recv`.

##### **Broadcasting Neighbor Maps**
- The `MPIHandler::broadcastBlockNeighborMap` method broadcasts the block adjacency map (neighbor relationships between blocks) from Rank 0 to all other ranks using `MPI_Bcast`.
- This ensures that all processes have the same view of the block relationships.

##### **Gathering Results**
- The `MPIHandler::gatherResults` method collects simulation results from all processes.
- Each process sends its local results (e.g., time step data for its cells) to Rank 0 using `MPI_Gatherv`.
- Rank 0 aggregates the results and normalizes the data (e.g., ensuring that the sum of S, I, and R values equals 1 for each cell).

---

#### 3️⃣ Rank Correlation with Code Implementation

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

#### 4️⃣ How Everything Is Put Together

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

#### 5️⃣ Scoring Logic for Dividing Data

The `GridSimulation::divideIntoOptimalBlocks` method evaluates multiple configurations for dividing the grid into blocks and selects the best one based on a scoring system. Here's how the scoring works:

##### **1. Configurations Evaluated**
- The method considers all divisors of the total number of cells as potential numbers of blocks.
- For each divisor (number of blocks), it calculates the number of cells per block.

##### **2. Factors Considered in Scoring**
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

##### **3. Scoring Formula**
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

##### **4. Example Output**
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

##### **5. Selected Configuration**
The configuration with the highest score is selected as the optimal distribution. In this example, the configuration with 10 blocks and 5 cells per block is chosen.

---

## Implementation Details

- The simulation uses a 2D grid of cells, where each cell represents a unit of the population.
- Each cell has SIR (Susceptible-Infected-Recovered) state variables, representing the number of susceptible, infected, and recovered individuals.
- The model updates the state of each cell at each time step based on the SIR dynamics and the interactions with neighboring cells.

---

## Contributors
- [Salvatore Mariano Librici](https://www.linkedin.com/in/salvatore-mariano-librici/)
- [Nada Khaled](https://www.linkedin.com/in/nadahkhaledd10/)
- [Milica Sanjevic](https://www.linkedin.com/in/milica-sanjevic-321392327/)
- [Hirdesh Kumar](https://www.linkedin.com/in/hirdeshkumar2407/)
- [Yibo Li](https://www.linkedin.com/in/yibo-li-0b0b792a2/)

## Analyzing Results

### Running Analysis
```bash
# After running both main simulation and tests
python scripts/analyze_results.py
```

### Output Location
Analysis results are saved in `./data/analysis/`:
- `sir_comparison.png`: Visual comparison of main simulation vs test results
- `error_metrics.csv`: Detailed error measurements
- `error_metrics.png`: Visualization of error metrics

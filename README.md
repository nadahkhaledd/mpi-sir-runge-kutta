# 🦠 PARALLEL SIR SIMULATION USING MPI & RUNGE-KUTTA METHOD.

A parallel implementation of the SIR (Susceptible-Infected-Recovered) model using MPI to simulate disease spread across The United States of America states.

## Data Source & Preprocessing

### Original Data Source
The datasets used in this project are sourced from the [JHU CSSE COVID-19 Dataset](https://github.com/CSSEGISandData/COVID-19/tree/4360e50239b4eb6b22f3a1759323748f36752177/csse_covid_19_data/csse_covid_19_daily_reports_us). Our main training dataset is from February 2, 2021.

### Data Preprocessing Steps
Before using a CSV file in the simulation:
1. Population data is added from external sources for each state
2. Missing values are filled using values from previous records
3. Dataset is cleaned to retain only the first 9 essential columns:
   - Province_State
   - Population
   - Last_Update
   - Lat
   - Long_
   - Confirmed
   - Deaths
   - Recovered
   - Active
4. States are reordered to place geographically adjacent states together
5. Header row is added with column count and row count

## Project Structure
```
.
├── data
│   ├── output                # Simulation results
│   ├── test_results           # Test outputs
│   ├── analysis              # Analysis plots and metrics
│   └── test_datasets         # Raw CSVs for testing
├── header
│    ├── main
│    └── test     
├── scripts                   # Python scripts for analysis and plotting
└── src                       # C++ source code
    ├── main.cpp              # Main simulation file
    ├── main
    └── test                  # Test suite for simulation
```

## Implementation Details

### Key Components
1. Data Distribution
   - Optimal block division
   - Load balancing
   - Neighbor cell mapping

2. MPI Communication
   - Block distribution
   - Ghost cell updates
   - Result gathering

3. SIR Model
   - Differential equations
   - Parameter tuning
   - State management

4. Output Handling
   - CSV writing
   - Logging
   - Error handling

### Performance Optimization
- Load balancing strategies
- Asynchronous communication
- Efficient I/O operations

## Building & Running

### Prerequisites
1. C++17 or later compiler
2. MPI Library
3. Python 3.x with required packages

### Installation Steps
```bash
# C++ Requirements
sudo apt-get update
sudo apt-get install build-essential
sudo apt-get install openmpi-bin libopenmpi-dev

# Python Requirements
sudo apt-get install python3 python3-pip
pip3 install numpy pandas matplotlib seaborn
```

### Main Simulation
```bash
# Create necessary directories
mkdir -p data/output

# Build
make clean
make all

# Run simulation (e.g., with 4 processes)
mpirun -np 4 ./sir_simulation

# Plot simulation results
python scripts/PlottingSIRModelResults.py
```

### Test Suite
```bash
# Create test directories
mkdir -p data/test_results data/analysis

# Build tests
make clean
make test

# Run tests
mpirun -np 4 ./sir_test_suite

# Analyze test results
python scripts/analyze_results.py
```

## Testing & Analysis

### Adding New Test Data
1. Place raw CSV in `data/test_datasets/`
2. Run preprocessing:
   ```bash
   python scripts/clean_sort_dataset.py data/test_datasets/your_dataset.csv
   ```
3. Verify the preprocessing steps:
   - Population data added
   - Missing values filled
   - Only essential columns retained
   - States geographically sorted
4. Update test configurations in `src/test/TestSuite.cpp`

### Available Test Datasets
1. sorted_01-01-2021.csv
   - First wave 2020 data
   - 50 states complete data
   - Used for base temporal tests

2. sorted_02-05-2021.csv
   - Second wave 2021 data
   - 50 states complete data
   - Used for comparative analysis

### Test Requirements
- Must contain all 50 US states
- Population values must be positive
- Missing values handled as zeros
- Dates in YYYY-MM-DD format

## Output & Analysis

### File Structure
- `data/output/`: Main simulation results
- `data/test_results/`: Test outputs
- `data/analysis/`: Analysis plots and metrics

### Analysis Scripts
1. Main Results:
   ```bash
   python scripts/PlottingSIRModelResults.py
   ```
2. Test Analysis:
   ```bash
   python scripts/analyze_results.py
   ```

## Simulation Results

### Output File Formats

#### Main Simulation Results
Location: `data/output/simulation_results.csv`
```csv
Time,S_avg,I_avg,R_avg
0.0,0.950000,0.050000,0.000000
0.2,0.947331,0.052669,0.000000
0.4,0.944516,0.055484,0.000000
...
```
Where:
- `Time`: Simulation timestep
- `S_avg`: Proportion of susceptible population
- `I_avg`: Proportion of infected population
- `R_avg`: Proportion of recovered population

#### Test Results
Location: `data/test_results/<test_name>_p<num_processes>_results.csv`
```csv
Time,S_avg,I_avg,R_avg
0.0,0.950000,0.050000,0.000000
...
```

#### Performance Metrics
Location: `data/output/timing_log.csv`
```csv
PhaseName,Statistic,Value,Units,NumRanks
distributeBlocks_Total,Min,0.000123,s,4
distributeBlocks_Total,Max,0.000145,s,4
distributeBlocks_Total,Avg,0.000134,s,4
...
```

### Generated Plots

#### 1. SIR Evolution
Location: `plots/sir_global_line_plot.png`
- Shows the temporal evolution of S, I, R populations
- X-axis: Time steps
- Y-axis: Population proportions
- Three lines: Susceptible (blue), Infected (red), Recovered (green)

#### 2. Infection Heatmap
Location: `plots/infection_heatmap_per_rank.png`
- Visualizes infection spread across MPI ranks
- X-axis: Time steps
- Y-axis: MPI ranks
- Color intensity: Infection level (darker = higher infection)

#### 3. Performance Analysis
Location: `plots/timing_comparison_phases.png`
- Compares execution times across simulation phases
- Shows min/max/avg times for each phase
- Helps identify performance bottlenecks

### Interpreting Results

1. **Convergence Check**
   - S + I + R should always sum to 1.0
   - Values should stabilize over time
   - Final R value indicates total affected population

2. **Performance Metrics**
   - Load balance: Compare execution times across ranks
   - Communication overhead: Check MPI phase timings
   - Scalability: Compare timings with different process counts

3. **Validation Criteria**
   - Infection peak timing matches historical data
   - Recovery rates align with medical observations
   - Geographic spread patterns are realistic

## 🔍 Detailed Implementation

### SIR Model Equations
The SIR model is based on the following set of differential equations:

<img src="https://latex.codecogs.com/svg.latex?\begin{align*}%20\frac{dS}{dt}%20&=%20-\beta%20\frac{SI}{N}%20\\%20\frac{dI}{dt}%20&=%20\beta%20\frac{SI}{N}%20-%20\gamma%20I%20\\%20\frac{dR}{dt}%20&=%20\gamma%20I%20\end{align*}" />

where:
- S, I, and R are the numbers of susceptible, infected, and recovered individuals
- N is the total population size (assumed constant)
- β (beta) is the transmission rate
- γ (gamma) is the recovery rate

### Parameter Tuning
Parameters are tuned based on:
- Literature values
- Calibration with observed data
- Sensitivity analysis to assess impact

### State Management
States are managed using a discrete event simulation approach:
- Events are scheduled for infections, recoveries, and data logging.
- Future events are predicted based on current state and parameters.
- State is updated at each event, and new events are scheduled as needed.

## Contributors
- [Salvatore Mariano Librici](https://www.linkedin.com/in/salvatore-mariano-librici/)
- [Nada Khaled](https://www.linkedin.com/in/nadahkhaledd10/)
- [Milica Sanjevic](https://www.linkedin.com/in/milica-sanjevic-321392327/)
- [Hirdesh Kumar](https://www.linkedin.com/in/hirdeshkumar2407/)
- [Yibo Li](https://www.linkedin.com/in/yibo-li-0b0b792a2/)

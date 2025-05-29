# Test Datasets Documentation

## Data Source
All datasets are sourced from [JHU CSSE COVID-19 Dataset](https://github.com/CSSEGISandData/COVID-19/tree/4360e50239b4eb6b22f3a1759323748f36752177/csse_covid_19_data/csse_covid_19_daily_reports_us)

## Preprocessing Steps
Each dataset undergoes the following preprocessing:
1. Population data is added from external sources for each state
2. Missing values are filled using values from previous records
3. States are reordered geographically for optimal simulation
4. Dataset is cleaned to retain only essential columns
5. Header row is added with dimensions

## Available Test Datasets

### Main Training Dataset
- **File**: sorted_02-02-2021.csv
- **Date**: February 2, 2021
- **Purpose**: Main simulation reference
- **Location**: `./data/sorted_initial_conditions.csv`

### Temporal Test Datasets
1. **Early Wave Dataset**
   - **File**: sorted_01-01-2021.csv
   - **Date**: January 1, 2021
   - **Purpose**: Testing model on early wave data
   - **Parameters**: β=0.3, γ=0.1

2. **Peak Wave Dataset**
   - **File**: sorted_02-05-2021.csv
   - **Date**: February 5, 2021
   - **Purpose**: Testing model on peak infection period
   - **Parameters**: β=0.3, γ=0.1

### Parameter Sensitivity Tests
All sensitivity tests use sorted_01-01-2021.csv with varying parameters:
1. **Low Transmission**
   - β=0.1
   - γ=0.1
   - Purpose: Test low infection spread

2. **High Transmission**
   - β=0.5
   - γ=0.1
   - Purpose: Test rapid infection spread

## File Format
Each CSV must contain:
```csv
num_rows num_cols
Province_State,Population,Last_Update,Lat,Long_,Confirmed,Deaths,Recovered,Active
Alabama,4903185,2021-01-01,...
...
```

## Required Columns
1. Province_State: US state name
2. Population: Total state population
3. Last_Update: Date in YYYY-MM-DD format
4. Lat: Latitude coordinate
5. Long_: Longitude coordinate
6. Confirmed: Total confirmed cases
7. Deaths: Total deaths
8. Recovered: Total recovered cases
9. Active: Currently active cases

## Data Validation
- Population > 0
- Confirmed ≥ Deaths + Recovered
- Active = Confirmed - Deaths - Recovered
- All 50 US states present
- No missing values (filled with previous values if missing)

## Adding New Test Data
1. Download data from JHU CSSE repository
2. Add population column with corresponding correct numbers
3. fill empty cells
4. Place raw CSV in this directory
5. Run preprocessing:
   ```bash
   python ../../scripts/clean_sort_dataset.py ./new_dataset.csv
   ```
6. Verify output format matches requirements
7. Add test configuration in `src/test/TestSuite.cpp`

## Directory Structure
```
data/
├── README.md
├── sorted_01-01-2021.csv
├── sorted_02-02-2021.csv
└── sorted_02-05-2021.csv
scripts/
└── clean_sort_dataset.py
src/
└── test/
    └── TestSuite.cpp
```

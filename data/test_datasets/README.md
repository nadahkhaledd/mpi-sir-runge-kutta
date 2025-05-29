# Test Datasets

## Dataset Format
Each dataset should include:
- First line: "<num_rows> <num_columns>"
- Second line: Column headers
- Data rows with columns:
  - Province_State
  - Population
  - Last_Update
  - Lat
  - Long_
  - Confirmed
  - Deaths
  - Recovered
  - Active

## Available Datasets
1. sorted_01-01-2021.csv
   - Period: First wave 2020
   - Size: 50 states
   - Source: CDC data
   - Used for: Base temporal tests

2. sorted_02-05-2021.csv
   - Period: Second wave 2021
   - Size: 50 states
   - Source: CDC data
   - Used for: Comparative analysis

## Processing New Datasets
1. Place raw CSV in this directory
2. Run preprocessing script:
   ```bash
   python ../../scripts/clean_sort_dataset.py ./raw_dataset.csv
   ```
3. Verify output format matches requirements
4. Update this README
5. Add test configuration in TestSuite.cpp

## Dataset Requirements
- Must contain all 50 US states
- Population values must be positive
- Missing values handled as zeros
- Dates in YYYY-MM-DD format

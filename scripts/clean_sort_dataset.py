import csv
import sys
import os

#python scripts/clean_sort_dataset.py ./data/initial_conditions.csv
# Define a geographical grid layout for US states
geographical_grid = [
    ["Washington", "Oregon", "Idaho", "Montana"],
    ["California", "Nevada", "Utah", "Wyoming"],
    ["Arizona", "New Mexico", "Colorado", "North Dakota"],
    ["Texas", "Oklahoma", "Kansas", "South Dakota"],
    ["Minnesota", "Iowa", "Missouri", "Nebraska"],
    ["Wisconsin", "Illinois", "Indiana", "Michigan"],
    ["Ohio", "Kentucky", "Tennessee", "West Virginia"],
    ["Pennsylvania", "New York", "Virginia", "North Carolina"],
    ["South Carolina", "Georgia", "Florida", "Alabama"],
    ["Arkansas", "Louisiana", "Mississippi", "Delaware"],
    ["Maryland", "New Jersey", "Connecticut", "Rhode Island"],
    ["Massachusetts", "Vermont", "New Hampshire", "Maine"]
]

# Flatten the grid into a single list of states in geographical order
state_order = [state for row in geographical_grid for state in row]

# List of valid US states to filter the CSV
us_states = [
    "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut",
    "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
    "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan",
    "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
    "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio",
    "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
    "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington", "West Virginia",
    "Wisconsin", "Wyoming"
]

def clean_sort_csv(input_file):
    input_dir = os.path.dirname(input_file)
    input_filename = os.path.basename(input_file)
    output_filename = f"sorted_{input_filename}"
    output_file = os.path.join(input_dir, output_filename)

    keep_columns = [0, 1, 3, 4, 5, 6, 7, 8, 9]
    column_names = [
        "Province_State", "Population", "Last_Update", 
        "Lat", "Long_", "Confirmed", "Deaths", "Recovered", "Active"
    ]

    # Read and process the CSV file
    with open(input_file, "r") as infile:
        reader = csv.reader(infile)
        header = next(reader)  # Skip original header
        rows = list(reader)

    # Filter rows to include only valid US states and keep only specified columns
    filtered_rows = []
    for row in rows:
        if row[0] in us_states:
            filtered_row = [row[i] for i in keep_columns]
            filtered_rows.append(filtered_row)

    # Create mapping of states to their data rows
    state_to_row = {row[0]: row for row in filtered_rows}

    # Sort rows based on geographical order
    sorted_rows = [state_to_row[state] for state in state_order if state in state_to_row]
    remaining_states = [row for row in filtered_rows if row[0] not in state_order]
    sorted_rows.extend(remaining_states)

    # Write the cleaned and sorted data
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        # Write dimensions line without comma (space-separated)
        outfile.write(f"{len(sorted_rows)} {len(keep_columns)}\n")
        writer.writerow(column_names)  # Write header
        writer.writerows(sorted_rows)  # Write data

    print(f"Cleaned and sorted CSV file saved to {output_file}")
    # print(f"Kept columns: {', '.join(column_names)}")
    # print(f"Total rows: {len(sorted_rows)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python clean_sort_dataset.py <input_csv_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.")
        sys.exit(1)
    
    clean_sort_csv(input_file)

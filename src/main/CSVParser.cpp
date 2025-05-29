#include "../../header/main/CSVParser.h"
#include "../../header/main/SIRCell.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include <unordered_map>
#include <map>
#include <algorithm>

std::string CSVParser::trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

std::vector<std::vector<double>> CSVParser::loadUSStateData(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error opening file " << filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    std::string line;
    bool headerFound = false;
    int lineCount = 0;
    
    // Each row will be stored as a vector of double
    std::vector<std::vector<double>> data;
    
    // Process lines until we find valid data
    while (std::getline(infile, line)) {
        lineCount++;
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Check if this line appears to be the header
        if (line.find("Province_State") != std::string::npos) {
            std::cout << "Found header at line " << lineCount << ": " << line << std::endl;
            headerFound = true;
            continue;
        }
        
        // Skip any lines before the header
        if (!headerFound) {
            std::cout << "Skipping pre-header line: " << line << std::endl;
            continue;
        }
        
        // Process data line
        std::istringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        // Split line on commas
        while (std::getline(ss, token, ',')) {
            tokens.push_back(trim(token));
        }
        
        // Verify we have enough columns
        if (tokens.size() < 9) {
            std::cerr << "Line " << lineCount << " has fewer than 9 columns: " << line << std::endl;
            continue;
        }
        
        try {
            // Convert values to doubles (not integers)
            double population = std::stod(tokens[1]); // Convert to double instead of int
            double lat = std::stod(tokens[3]);
            double lon = std::stod(tokens[4]);
            double confirmed = std::stod(tokens[5]);
            double deaths = std::stod(tokens[6]);
            
            // Handle empty or missing Recovered/Active fields
            double recovered = tokens[7].empty() ? 0.0 : std::stod(tokens[7]);
            double active = tokens[8].empty() ? confirmed - deaths - recovered : std::stod(tokens[8]);
            
            // Store values
            data.push_back({population, lat, lon, confirmed, deaths, recovered, active});
        } catch (const std::exception& e) {
            std::cerr << "Invalid value at line " << lineCount << ": " << line << "\nError: " << e.what() << std::endl;
            continue;
        }
    }
    
    std::cout << "Successfully parsed " << data.size() << " data rows from CSV." << std::endl;
    return data;
}

SIRCell CSVParser::mapToSIR(const std::vector<double>& rowData) {
    // Adjusted to match the correct column structure:
    // Province_State, Population, Date, Lat, Long, Confirmed, Deaths, Recovered, Active
    double population = rowData[0]; // Population
    double confirmed = rowData[3]; // Confirmed cases
    double deaths = rowData[4];    // Deaths
    double recovered = rowData[5]; // Recovered
    double active = rowData[6];    // Active cases

    // Ensure population is valid
    if (population <= 0) {
        std::cerr << "Warning: Invalid population value, using default values\n";
        return SIRCell(1.0, 0.0, 0.0);  // Return default state
    }

    // Calculate S, I, R with safety checks
    double S = (population - confirmed) / population;
    double I = (active > 0) ? active / population : 0.0;
    double R = ((recovered + deaths) > 0) ? (recovered + deaths) / population : 0.0;

    // Normalize to ensure sum is 1.0
    double sum = S + I + R;
    if (sum > 0) {
        S /= sum;
        I /= sum;
        R /= sum;
    } else {
        S = 1.0;
        I = 0.0;
        R = 0.0;
    }

    return SIRCell(S, I, R);
}
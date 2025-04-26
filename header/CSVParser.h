#ifndef CSVPARSER_H
#define CSVPARSER_H

#include <string>
#include <vector>
#include "SIRCell.h"

class CSVParser {
private:
    // Helper function: trim whitespace from a string
    static std::string trim(const std::string &s);

public:
    // Parse CSV data
    static std::vector<std::vector<double>> loadUSStateData(const std::string& filename);
};

#endif // CSVPARSER_H
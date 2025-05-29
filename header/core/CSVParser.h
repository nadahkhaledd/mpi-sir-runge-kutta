#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include "SIRCell.h"
#include <string>
#include <vector>
#include <stdexcept>

class CSVParser {
public:
    static std::string trim(const std::string& s);
    static std::vector<std::vector<double>> loadUSStateData(const std::string& filename);
    static SIRCell mapToSIR(const std::vector<double>& rowData);
    
private:
    static double safeStod(const std::string& str, double defaultValue = 0.0);
};

#endif
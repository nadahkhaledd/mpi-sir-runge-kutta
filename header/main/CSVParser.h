#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include "SIRCell.h"
#include <string>
#include <vector>

class CSVParser {
public:
    static std::string trim(const std::string& s);
    static std::vector<std::vector<double>> loadUSStateData(const std::string& filename);
    static SIRCell mapToSIR(const std::vector<double>& rowData);
};

#endif
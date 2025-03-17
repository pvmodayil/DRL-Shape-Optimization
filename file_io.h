#ifndef FILE_IO_H
#define FILE_IO_H

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

namespace fileio{

// Function to split a string into a vector of strings based on a delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

// Function to read a CSV file and store each column as a vector
std::unordered_map<std::string, std::vector<double>> readCSV(const std::string& filename);

} // namespace fileio

#endif
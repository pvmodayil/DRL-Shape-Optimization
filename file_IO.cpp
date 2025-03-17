#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>

// Function to split a string into a vector of strings based on a delimiter
std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

// Function to read a CSV file and store each column as a vector
std::unordered_map<std::string, std::vector<double>> readCSV(const std::string& filename) {
    // Create map to be returned
    std::unordered_map<std::string, std::vector<double>> data;
    // Check if the file can be opened
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file " << filename << std::endl;
        return;
    }

    // Initialize the line
    std::string line;
    // Split the header line into individual header_names
    std::getline(file, line);
    std::vector<std::string> column_headers = split(line, ',');

    // Add the column names and create empty vectors
    for (size_t i = 0; i < column_headers.size(); ++i) {
        data.insert({column_headers[i],std::vector<double>()});
    }

    // Read rest of the lines in the file
    while (std::getline(file, line)) {
        std::vector<std::string> row = split(line, ',');

        // Check if the row has the same number of header_names as the header
        if (row.size() != column_headers.size()) {
            std::cerr << "Error: Row has " << row.size() << "values, for " << column_headers.size() << " headers" << std::endl;
            continue;
        }
        // Store each column as a vector
        for (size_t i = 0; i < column_headers.size(); ++i) {
            data[column_headers[i]].push_back(std::stod(row[i])); // convert the string into double values
        }
    }

    file.close();

    return data;
}
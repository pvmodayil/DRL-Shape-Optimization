#include <file_io.h>
#include <iostream>
#include <vector>
#include "microstrip_arrangement.h"

int main(){
    
    std::string filename = "result_curve.csv";
    
    std::unordered_map<std::string, std::vector<double>> data = fileio::readCSV(filename);
    // Print the data to verify
    for (const auto& pair : data) {
        std::cout << pair.first << ": ";
        for (const auto& value : pair.second) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
#include <file_io.h>
#include <microstrip_arrangement.h>
#include <iostream>
#include <vector>
#include "microstrip_arrangement.h"

struct MicrostripArrangement {
    const double V0;
    const double hw_micrstr;
    const double ht_micrstr;
    const double hw_arra;
    const double ht_arra;
    const double ht_subs;
    const double er1;
    const double er2;
    const int N; 
};
struct GeneticAlgo {
    const int N;
    const int evolution_steps;
    const int population_size;
    const int noise_scale;
};

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

    MicrostripArrangement arrangement = {
        .V0 = 1.0,
        .hw_micrstr = 0.05e-3,  // Half-width of microstrip in meters 
        .ht_micrstr = 0.0,   // Height of microstrip in meters 
        .hw_arra = 1.38e-3,      // Half-width of array in meters 
        .ht_arra = 2.76e-3,       // Height of array in meters 
        .ht_subs = 0.1382e-3,     // Height of substrate in meters 
        .er1 = 1.0,            // Relative permittivity of substrate
        .er2 = 12.9,            // Relative permittivity of air
        .N = 2000                 // Number of elements in the array
    };

    // Filter vectors
    MSA::filterVectors(arrangement.hw_micrstr,
                    arrangement.hw_arra,
                    data["g_ptsy"],
                    data["g_ptsx"],
                    data["gpts_y"],
                    data["g_ptsx"]);
                    
    std::cout << "Checking for monotonicity" << std::endl;
    std::cout << MSA::isMonotonicallyDecreasing(data["g_ptsy"]) << std::endl;
    
    return 0;
}
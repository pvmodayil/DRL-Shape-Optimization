#include <file_io.h>
#include <microstrip_arrangement.h>
#include<genetic_algorithm.h>
#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>

// External Includes
#include <Eigen/Dense>

int main(){
    std::cout << "Genetic Algorithm Optimization" << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads." << std::endl;
    // Set the number of threads for OpenMP
    omp_set_num_threads(omp_get_max_threads());
    // Read starting curve
    std::string filename = "../result_curve.csv";
    std::cout<< "Reading g point values from: " << filename << std::endl;
    std::unordered_map<std::string, std::vector<double>> data = fileio::readCSV(filename);

    // Microstrip arrangement 
    MSA::MicrostripArrangement arrangement = MSA::MicrostripArrangement(
        1.0, // V0
        0.05e-3, // d
        0.0, // t
        1.38e-3, // a
        2.76e-3, // b
        0.1382e-3,// c/h
        1.0, // er1
        12.9, // er2
        2000); // N

    std::vector<double> x = data["g_ptsx"];
    std::vector<double> g = data["g_ptsy"];
    // Filter the vectors if necessary
    if(x[0] <= arrangement.hw_micrstr || x.back() >= arrangement.hw_arra){
        // first pair of g,x are passed as value and the second is passed as reference so the original vectors itself will be filtered
        MSA::filterVectors(arrangement.hw_micrstr,arrangement.hw_arra,g,x,g,x);
    }

    // Convert the x and g vectors to Eigen arrays
    Eigen::ArrayXd x_array = Eigen::Map<const Eigen::ArrayXd>(x.data(), x.size(), 1); // Mx1
    Eigen::ArrayXd g_array = Eigen::Map<const Eigen::ArrayXd>(g.data(), g.size(), 1); // Mx1
    
    // Initial energy calculation
    Eigen::ArrayXd vn = MSA::calculatePotentialCoeffs(arrangement.V0,
        arrangement.hw_micrstr,
        arrangement.hw_arra,
        arrangement.N,
        g_array,
        x_array);

    double energy = MSA::calculateEnergy(arrangement.er1,
        arrangement.er2,
        arrangement.hw_arra,
        arrangement.ht_arra,
        arrangement.ht_subs,
        arrangement.N,
        vn);

    // Genetic Algorithm class
    int population_size = 100; // Number of individuals in the population
    int num_generations = 1000; // Number of generations to evolve
    GA::GeneticAlgorithm ga_problem = GA::GeneticAlgorithm(arrangement,g_array,x_array,population_size,num_generations,0.1);
    double noise_scale = 0.1;

    // Result struct
    GA::GeneticResult result = GA::GeneticResult(Eigen::ArrayXd(num_generations+1), Eigen::VectorXd(g.size() + 3), 0.0); // Empty initialization
    result.energy_convergence(0) = energy; // Initial energy
    
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    // Call the function you want to time
    ga_problem.optimize(noise_scale, result);
    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();
    
    std::cout << "Best energy: " << result.best_energy << "\n";
    std::cout << "Best curve: " << result.best_curve << "\n";

    // Calculate the duration
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Execution time: " << duration.count()/1000 << " s" << std::endl;

    // Convert the result to a vector and save it to a CSV file
    std::vector<double> result_vector(result.best_curve.data(), result.best_curve.data() + result.best_curve.size());
    std::unordered_map<std::string, std::vector<double>> result_data;
    result_data["g_ptsx"] = data["g_ptsx"];
    result_data["g_ptsy"] = result_vector;
    std::string curve_output_filename = "../result_curve_optimized.csv";
    fileio::writeCSV(curve_output_filename, result_data);

    std::unordered_map<std::string, std::vector<double>> energy_history;
    energy_history["energy"] = std::vector<double>(result.energy_convergence.data(), result.energy_convergence.data() + result.energy_convergence.size());
    energy_history["generation"] = std::vector<double>(num_generations + 1);
    for (size_t i = 0; i < num_generations + 1; ++i) {
        energy_history["generation"][i] = static_cast<double>(i);
    }

    std::string history_output_filename = "../energy_history.csv";
    fileio::writeCSV(history_output_filename, energy_history);

    std::cout << "Energy history saved to: " << history_output_filename << std::endl;
    std::cout << "Optimized curve saved to: " << curve_output_filename << std::endl;
    return 0;
}
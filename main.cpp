#include <file_io.h>
#include <microstrip_arrangement.h>
#include<genetic_algorithm.h>
#include <iostream>
#include <vector>
#include <chrono>

// External Includes
#include <Eigen/Dense>

int main(){
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
    
    // Genetic Algorithm class
    GA::GeneticAlgorithm ga_problem = GA::GeneticAlgorithm(arrangement,g_array,x_array,100,1000,0.1);
    double noise_scale = 0.1;
    
    // Start timing
    auto start = std::chrono::high_resolution_clock::now();
    // Call the function you want to time
    ga_problem.optimize(noise_scale);
    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double, std::milli> duration = end - start;

    std::cout << "Execution time: " << duration.count()/1000 << " s" << std::endl;
    

    // Energy calculation
    Eigen::ArrayXd vn = MSA::calculatePotentialCoeffs(arrangement.V0,
        arrangement.hw_micrstr,
        arrangement.hw_arra,
        arrangement.N,
        g_array,
        x_array);
    
    Eigen::ArrayXd VF = MSA::calculatePotential(arrangement.hw_arra,
        arrangement.N,
        vn,
        data["g_ptsx"]);

    double energy = MSA::calculateEnergy(arrangement.er1,
        arrangement.er2,
        arrangement.hw_arra,
        arrangement.ht_arra,
        arrangement.ht_subs,
        arrangement.N,
        vn);
    
    std::cout<<"Energy: "<<energy<<std::endl;
    
    return 0;
}
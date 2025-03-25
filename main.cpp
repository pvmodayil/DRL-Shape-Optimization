#include <file_io.h>
#include <microstrip_arrangement.h>
#include<genetic_algorithm.h>
#include <iostream>
#include <vector>
struct GeneticAlgo {
    const int N;
    const int evolution_steps;
    const int population_size;
    const int noise_scale;
};

int main(){
    std::string filename = "../nominal.csv";
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
        5000); // N

    // Genetic Algorithm class
    GA::GeneticAlgorithm ga_problem = GA::GeneticAlgorithm(arrangement,data["g_ptsy"],data["g_ptsx"],100,0,0.1);
    double noise_scale = 0.1;
    ga_problem.optimize(noise_scale);

    // Energy calculation
    Eigen::ArrayXd vn = MSA::calculatePotentialCoeffs(arrangement.V0,
        arrangement.hw_micrstr,
        arrangement.hw_arra,
        arrangement.N,
        data["g_ptsy"],
        data["g_ptsx"]);
    
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
#include <file_io.h>
#include <microstrip_arrangement.h>
#include <iostream>
#include <vector>

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

    MicrostripArrangement(
        const double V0, 
        const double hw_micrstr, // Half-width of microstrip in meters 
        const double ht_micrstr, // Height of microstrip in meters 
        const double hw_arra, // Half-width of array in meters 
        const double ht_arra, // Height of array in meters 
        const double ht_subs, // Height of substrate in meters 
        const double er1, // Relative permittivity of air
        const double er2, // Relative permittivity of substrate
        const int N) // Number of Fourier coefficients
        : 
        V0(V0), 
        hw_micrstr(hw_micrstr), 
        ht_micrstr(ht_micrstr), 
        hw_arra(hw_arra), 
        ht_arra(ht_arra), 
        ht_subs(ht_subs), 
        er1(er1), 
        er2(er2), 
        N(N){}
};
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
    MicrostripArrangement arrangement = MicrostripArrangement(1.0,
        0.05e-3,
        0.0,
        1.38e-3,
        2.76e-3,
        0.1382e-3,
        1.0,
        12.9,
        5000);
        
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
    std::cout<<"Potential: "<<VF<<std::endl;

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
#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <Eigen/Dense>
#include <vector>
#include <random>

namespace GA{

    // Constructor
    GeneticAlgorithm::GeneticAlgorithm(const MSA::MicrostripArrangement& arrangement,
        std::vector<double>& starting_curveY,
        std::vector<double>& starting_curveX, 
        int population_size, 
        int num_generations, 
        double mutation_rate)
        : 
        arrangement(arrangement),
        starting_curveY(starting_curveY),
        starting_curveX(starting_curveX), 
        population_size(population_size), 
        num_generations(num_generations), 
        mutation_rate(mutation_rate) {}
    
    // Initialize the population with random noise
    Eigen::MatrixXd GeneticAlgorithm::initializePopulation(const double& noise_scale) const{
        const int vector_size = starting_curveY.size();
        const Eigen::ArrayXd column = Eigen::Map<const Eigen::ArrayXd>(starting_curveY.data(), vector_size, 1);
        
        // Create a zero matrix as init population
        Eigen::MatrixXd initialPopulation = Eigen::MatrixXd(vector_size, population_size);
        for (size_t i = 0; i < population_size; ++i) {
            initialPopulation.col(i) = column;
        }

        initialPopulation = initialPopulation + noise_scale * Eigen::MatrixXd::Random(vector_size, population_size);
        
        initialPopulation = initialPopulation.array().max(arrangement.V0).matrix();

        return initialPopulation;
    }

    double GeneticAlgorithm::calculateFitness(const std::vector<double>& individual) const{
        // Energy calculation
        Eigen::ArrayXd vn = MSA::calculatePotentialCoeffs(arrangement.V0,
            arrangement.hw_micrstr,
            arrangement.hw_arra,
            arrangement.N,
            individual,
            starting_curveX);

        return MSA::calculateEnergy(arrangement.er1,
            arrangement.er2,
            arrangement.hw_arra,
            arrangement.ht_arra,
            arrangement.ht_subs,
            arrangement.N,
            vn);
    }

    // Main function to run the optimization
    void GeneticAlgorithm::optimize(const double& noise_scale){
        // Create an initial population
        Eigen::MatrixXd initialPopulation = initializePopulation(noise_scale);
    }
    

} // end of namespace
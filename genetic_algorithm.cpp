#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <iostream>

namespace GA{

    // Constructor
    GeneticAlgorithm::GeneticAlgorithm(MSA::MicrostripArrangement& arrangement,
        Eigen::ArrayXd& starting_curveY,
        Eigen::ArrayXd& starting_curveX, 
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
        // Need the length for further processing
        int vector_size = starting_curveY.size();
        
        // Map the vector and create a matrix where each column is a copy of the starting curve
        Eigen::ArrayXd column = Eigen::Map<const Eigen::ArrayXd>(starting_curveY.data(), vector_size, 1);
        Eigen::MatrixXd initialPopulation = column.replicate(1, population_size);
        
        // Random uniform distribution between -1 to 1
        Eigen::MatrixXd random_matrix = Eigen::MatrixXd::Random(vector_size, population_size); // Think about scaling to 0 to 1

        Eigen::MatrixXd random_noise = (
            (noise_scale * random_matrix).array() * initialPopulation.array()
        ).matrix(); // do element wise multiplication to scale the noise for the starting curve

        // Add noise to create the initial population
        initialPopulation = initialPopulation + random_noise;

        // Limit the initial population within the boundary(V0)
        initialPopulation = initialPopulation.array().min(arrangement.V0).matrix();

        return initialPopulation;
    }

    double GeneticAlgorithm::calculateFitness(Eigen::ArrayXd& individual) const{
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

        // Fitness calculation

        // Random Crossover and Mutate

        //Repeat
        
    }
    

} // end of namespace
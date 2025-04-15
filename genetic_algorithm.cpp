#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric> // Include this header for std::iota
#include <iostream>
#include <stdexcept>
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
    Eigen::MatrixXd GeneticAlgorithm::initializePopulation(double& noise_scale){
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

    double GeneticAlgorithm::calculateFitness(Eigen::ArrayXd& individual){
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
    
    // Select the best and worst performers for Elitism implementation
    std::vector<size_t> GeneticAlgorithm::selectElits(const Eigen::ArrayXd& fitness_array){
        size_t n = fitness_array.size();
        if (n < 4){
            throw std::invalid_argument("Not enough individuals in population get at least 6!"); // Not enough elements
        }

        // Initialise the indexes
        size_t min1 = 0, min2 = 1, max1 = 0, max2 = 1;
    
        // Initialize min1, min2, max1, max2
        if (fitness_array[1] < fitness_array[0]) {
            std::swap(min1, min2);
        } else {
            std::swap(max1, max2);
        }
    
        for (size_t i = 2; i < n; ++i) {
            // Update for minimum values
            if (fitness_array[i] < fitness_array[min1]) {
                min2 = min1;
                min1 = i;
            } else if (fitness_array[i] < fitness_array[min2]) {
                min2 = i;
            }
    
            // Update for maximum values
            if (fitness_array[i] > fitness_array[max1]) {
                max2 = max1;
                max1 = i;
            } else if (fitness_array[i] > fitness_array[max2]) {
                max2 = i;
            }
        }
    
        return {min1, min2, max1, max2};
    }
    // Select the best and worst performers
    std::vector<size_t> GeneticAlgorithm::selectParents(const Eigen::ArrayXd& fitness_array) {
        // Get the best and worst performers
        std::vector<size_t> elits_indices = selectElits(fitness_array);

    }

    // Crossover
    Eigen::MatrixXd GeneticAlgorithm::crossover(Eigen::ArrayXd& parent1, Eigen::ArrayXd& parent2){
        // pass
        Eigen::MatrixXd empty;
        return empty;
    }

    //Mutation
    void GeneticAlgorithm::mutate(Eigen::ArrayXd& individual, double& noise_scale){
        // pass
    }

    // Main function to run the optimization
    void GeneticAlgorithm::optimize(double& noise_scale){
        
        // Create an initial population
        Eigen::MatrixXd population = initializePopulation(noise_scale); 
        
        // Iterate for num_generations steps
        for(size_t genration=0; genration<num_generations; ++genration){
            // Fitness calculation
            Eigen::ArrayXd fitness_array = Eigen::ArrayXd(population_size);
            for(size_t i =0; i<population_size; ++i){
                Eigen::ArrayXd individual = population.col(i);
                fitness_array[i] = calculateFitness(individual);
            }

            // Select potential parents and worst performers
            std::vector<size_t> selected_indices = selectParents(fitness_array);

            std::cout << "Fitness Array:\n" << fitness_array << std::endl;
            std::cout << " Top 2 Least values are: " << fitness_array[selected_indices[0]] << " , " << fitness_array[selected_indices[1]] << std::endl;
            std::cout << " Top 2 Largest values are: " << fitness_array[selected_indices[2]] << " , " << fitness_array[selected_indices[3]] << std::endl;
            Eigen::ArrayXd parent1 = population.col(selected_indices[0]);
            Eigen::ArrayXd parent2 = population.col(selected_indices[1]);

            //std::cout << "Fitness values: " << fitness_array << std::endl;
            // Random Crossover and Mutate

            //Repeat
            population = initializePopulation(noise_scale); // replace least performing 2 with children
        }
        
    }
    

} // end of namespace
#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric> // Include this header for std::iota
#include <iostream>
#include <stdexcept>
#include <omp.h>

void printProgressBar(int total, int current) {
    const int bar_width = 20; // Width of the progress bar
    float progress = static_cast<float>(current) / total;

    std::cout << "[";
    int pos = static_cast<int>(bar_width * progress);
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << static_cast<int>(progress * 100.0f) << "%\r"; // \r returns cursor to the beginning of the line
    std::cout.flush(); // Ensure the output is printed immediately
}

namespace GA{

    // Constructor
    // ------------------------------------------------------
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
    // ------------------------------------------------------
    Eigen::MatrixXd GeneticAlgorithm::initializePopulation(double& noise_scale){
        // Need the length for further processing
        int vector_size = starting_curveY.size();
        
        // Map the vector and create a matrix where each column is a copy of the starting curve
        Eigen::ArrayXd column = Eigen::Map<const Eigen::ArrayXd>(starting_curveY.data(), vector_size, 1);
        Eigen::MatrixXd initial_population = column.replicate(1, population_size);
        
        // Random uniform distribution between -1 to 1 scaled to 0 to 1 and further scaled to create random noise
        Eigen::MatrixXd random_noise = (
            (noise_scale * 
                0.5 * (Eigen::MatrixXd::Ones(vector_size, population_size) + Eigen::MatrixXd::Random(vector_size, population_size))
            ).array() * initial_population.array()
        ).matrix();

        // Add noise to create the initial population
        initial_population = initial_population + random_noise;

        // Limit the initial population within the boundary(V0)
        initial_population = initial_population.array().min(arrangement.V0).matrix();

        return initial_population;
    }

    // Fitness operator
    // ------------------------------------------------------
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
    
    // Selection operator
    // ------------------------------------------------------
    // Select the best and worst performers for Elitism implementation
    std::vector<size_t> GeneticAlgorithm::selectElites(const Eigen::ArrayXd& fitness_array){
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

    // Select the parents
    std::vector<size_t> GeneticAlgorithm::selectParents(const std::vector<size_t>& elites_indices, const Eigen::ArrayXd& fitness_array) {
        // Tournament selection
        const int TOURNAMENT_SIZE = 3;
        std::vector<size_t> selected_indices;
        std::vector<size_t> candidate_indices;
        size_t candidate_index;
        std::random_device rd; // random number from machine to put random seed
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, fitness_array.size() - 1);

        // Select top two out of three random selections
        candidate_indices.clear(); // Start fresh
        for (int i = 0; i < TOURNAMENT_SIZE; ++i){
            candidate_index = dis(gen);
            // // Do not require the worst performers as candidates and not 
            // while (candidate_index == elites_indices[2] || candidate_index == elites_indices[3]) {
            //     candidate_index = dis(gen);
            // }

            candidate_indices.push_back(candidate_index);
        }
        
        // Find the top two candidates out of the random three      
        size_t min1 = candidate_indices[0], min2 = candidate_indices[1];
        if (fitness_array[min2] < fitness_array[min1]) {
            std::swap(min1, min2);
        }
        if (fitness_array[candidate_indices[2]] < fitness_array[min1]) {
            min2 = min1;
            min1 = candidate_indices[2];
        } else if (fitness_array[candidate_indices[2]] < fitness_array[min2]) {
            min2 = candidate_indices[2];
        }
        
        selected_indices.push_back(min1);
        selected_indices.push_back(min2);

        return selected_indices;
    }

    // Crossover
    Eigen::VectorXd GeneticAlgorithm::crossover(Eigen::VectorXd& parent1, Eigen::VectorXd& parent2){
        
        size_t parent_size = parent1.size();
        // Random distribution initialize
        std::random_device rd; // random number from machine to put random seed
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, parent_size - 1);
        
        // Choose the random crossover point
        size_t crossover_point = dis(gen);

        // Crossover
        Eigen::VectorXd child = parent1;
        child.segment(crossover_point, parent_size - crossover_point) = parent2.segment(crossover_point, parent_size - crossover_point);

        // std::cout << "Parent 1:" << parent1 <<"\n";
        // std::cout << "Parent 2:" << parent2 <<"\n";
        // std::cout << "Child:" << child <<"\n";
        
        return child;
    }

    // Reproduction operator
    // ------------------------------------------------------
    Eigen::MatrixXd GeneticAlgorithm::reproduce(Eigen::MatrixXd& population, Eigen::ArrayXd& fitness_array, double& noise_scale, std::vector<double>& energies){
        // Inits
        std::vector<size_t> selected_indices;
        Eigen::VectorXd parent1;
        Eigen::VectorXd parent2;
        size_t vector_size = starting_curveY.size();

        // Create a random noise scaled matrix for mutation
        Eigen::MatrixXd new_population = ((Eigen::MatrixXd::Random(vector_size, population_size)*noise_scale).array() * 
        population.array()).matrix(); // do element wise multiplication to scale the noise for the curve

        // Get the best and worst performers
        std::vector<size_t> elites_indices = selectElites(fitness_array);
        energies.push_back(fitness_array[elites_indices[0]]);
        // std::cout << "Best Curve: " << population.col(elites_indices[0]) << "\n";
        
        // Reproduction cycle for population_size - 2 , need to retain the two elites
        // #pragma omp parallel for
        for (size_t i=0; i<population_size; ++i){
            selected_indices = selectParents(elites_indices, fitness_array);
            parent1 = population.col(selected_indices[0]);
            parent2 = population.col(selected_indices[1]);
            new_population.col(i) += crossover(parent1,parent2); // Crossover + mutate
        }
        
        // Retain the elites
        // new_population.col(population_size-2) = population.col(elites_indices[0]);
        // new_population.col(population_size-1) = population.col(elites_indices[1]);

        // Limit within the thresholds
        new_population = new_population.array().min(arrangement.V0).matrix();
        new_population = new_population.array().max(0).matrix();
        return new_population;
    }

    // Main function to run the optimization
    // ------------------------------------------------------
    void GeneticAlgorithm::optimize(double& noise_scale){
        
        std::vector<double> energies;
        // Create an initial population
        Eigen::MatrixXd population = initializePopulation(noise_scale); 
        Eigen::ArrayXd fitness_array = Eigen::ArrayXd(population_size);
        // Iterate for num_generations steps
        for(size_t genration=0; genration<num_generations; ++genration){
            printProgressBar(num_generations, genration+1);
            // Fitness calculation
            // #pragma omp parallel for
            for(size_t i =0; i<population_size; ++i){
                Eigen::ArrayXd individual = population.col(i);
                fitness_array[i] = calculateFitness(individual);
            }

            // Reproduce
            population = reproduce(population, fitness_array, noise_scale, energies);
        }
        std::cout << "Energy values: ";
        for (auto energy:energies){
            std::cout << energy << "," ;
        }
        std::cout << "\n";

        std::vector<size_t> elites_indices = selectElites(fitness_array);
        std::cout << "\nBest Energy: " << fitness_array[elites_indices[0]] << "\n";
        std::cout << "Best Curve: " << population.col(elites_indices[0]) << "\n";
    }
    

} // end of namespace
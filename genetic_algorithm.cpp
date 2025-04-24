#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <tuple> // For std::tuple
#include <iostream>
#include <stdexcept>

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
        mutation_rate(mutation_rate) {
            std::random_device rd;
            gen = std::mt19937(rd());
            real_dist = std::uniform_real_distribution<>(0.0, 1.0);
            parent_index_dist = std::uniform_int_distribution<>(0, population_size - 1);
        }
    
    // Initialize the population with random noise
    // ------------------------------------------------------
    Eigen::MatrixXd GeneticAlgorithm::initializePopulation(double& noise_scale){
        // Need the length for further processing
        size_t vector_size = starting_curveY.size();
        
        // Map the vector and create a matrix where each column is a copy of the starting curve
        Eigen::MatrixXd initial_population = starting_curveY.replicate(1, population_size);
        
        // Random uniform distribution between -1 to 1 scaled to 0 to 1 and further scaled to create random noise
        Eigen::MatrixXd random_noise = (
            (noise_scale * 
                0.5 * (Eigen::MatrixXd::Ones(vector_size, population_size) + Eigen::MatrixXd::Random(vector_size, population_size))
            ).array() * initial_population.array()
        ).matrix();

        // Add noise to create the initial population
        initial_population.noalias() += random_noise;

        // Limit the initial population within the boundary(V0)
        initial_population = initial_population.array().min(arrangement.V0).max(0).matrix();

        return initial_population;
    }

    // Fitness operator
    // ------------------------------------------------------
    double GeneticAlgorithm::calculateFitness(Eigen::ArrayXd& individual){
        if (!MSA::isMonotonicallyDecreasing(individual)){
            return 100.0; // high energy value since necessary condition failed
        }
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
    size_t GeneticAlgorithm::selectElites(const Eigen::ArrayXd& fitness_array){
        size_t n = fitness_array.size();

        // Initialise the indexes
        size_t min = 0;
    
        for (size_t i = 0; i < n; ++i) {
            // Update for minimum values
            if (fitness_array[i] < fitness_array[min]) {
                min = i;
            }
        }
    
        return min;
    }

    // Select the parents
    size_t GeneticAlgorithm::selectParent(const Eigen::ArrayXd& fitness_array) {
        // Tournament selection with size 2
        size_t candidate_index1;
        size_t candidate_index2;

        // Do tournament selection
        candidate_index1 = parent_index_dist(gen);
        candidate_index2 = parent_index_dist(gen);

        if (fitness_array[candidate_index2] < fitness_array[candidate_index1]) {
            std::swap(candidate_index1, candidate_index2);
        }

        return candidate_index1;
    }

    // Crossover
    std::tuple<Eigen::VectorXd, Eigen::VectorXd> GeneticAlgorithm::crossover(Eigen::VectorXd& parent1, Eigen::VectorXd& parent2, double eta){
        
        size_t parent_size = parent1.size();
        double exponent = 1.0 / (eta + 1.0);

        // Generate a vector of random numbers
        //Eigen::ArrayXd u = Eigen::ArrayXd::NullaryExpr(parent_size, [this]() { return real_dist(gen); }); // pass the class instance
        Eigen::ArrayXd u = 0.5 * (Eigen::ArrayXd::Random(parent_size) + 1); // Random numbers between 0 and 1

        // The following arithmetic has 1.0 - u in the denominator and hence it is important that u never has 1.0 as value
        double epsilon = std::numeric_limits<double>::epsilon();
        u = u.cwiseMax(epsilon).cwiseMin(1.0 - epsilon); // Ensure u is not 1.0 or 0.0

        Eigen::ArrayXd beta(u.size());
        Eigen::Array<bool, Eigen::Dynamic, 1> mask = (u <= 0.5); // mask array
        // Compute both cases for the mask
        Eigen::ArrayXd beta_case_ulessthan  = (2.0 * u).pow(exponent);
        Eigen::ArrayXd beta_case_umorethan = (1.0 / (2.0 * (1.0 - u))).pow(exponent);
        beta = mask.select(beta_case_ulessthan, beta_case_umorethan);

        // Calculate children
        Eigen::VectorXd child1 = 0.5 * ((1.0 + beta) * parent1.array() + (1.0 - beta) * parent2.array());
        Eigen::VectorXd child2 = 0.5 * ((1.0 - beta) * parent1.array() + (1.0 + beta) * parent2.array());

        return std::make_tuple(child1, child2);
    }

    // Reproduction operator
    // ------------------------------------------------------
    Eigen::MatrixXd GeneticAlgorithm::reproduce(Eigen::MatrixXd& population, Eigen::ArrayXd& fitness_array, double& noise_scale){
        // Inits
        std::vector<size_t> selected_indices;
        size_t vector_size = starting_curveY.size();
        Eigen::VectorXd parent1(vector_size);
        Eigen::VectorXd parent2(vector_size);

        // Create a random noise scaled matrix for mutation
        Eigen::MatrixXd new_population(vector_size, population_size);
        
        for (size_t i=0; i<population_size; i+=2){
            parent1 = population.col(selectParent(fitness_array));
            parent2 = population.col(selectParent(fitness_array));
            auto [child1, child2] = crossover(parent1,parent2); // Crossover + mutate
            new_population.col(i).noalias() = child1;
            new_population.col(i+1).noalias() = child2;

        }

        // Limit within the thresholds
        new_population = new_population.array().min(arrangement.V0).matrix();
        new_population = new_population.array().max(0).matrix();
        return new_population;
    }

    // Main function to run the optimization
    // ------------------------------------------------------
    void GeneticAlgorithm::optimize(double& noise_scale){
        
        size_t vector_size = starting_curveY.size();
        // Create an initial population
        Eigen::MatrixXd population = initializePopulation(noise_scale); 
        Eigen::ArrayXd fitness_array = Eigen::ArrayXd(population_size);
        // Eigen::ArrayXd individual;
        // Iterate for num_generations steps
        for(size_t generation=0; generation<num_generations; ++generation){
            printProgressBar(num_generations, generation + 1);

            // Fitness calculation
            #pragma omp parallel for
            for(int i=0; i<population_size; ++i){
                Eigen::ArrayXd individual = population.col(i).array();
                fitness_array[i] = calculateFitness(individual);
            }

            // Reproduce
            population = reproduce(population, fitness_array, noise_scale);
        }

        size_t elites_index = selectElites(fitness_array);
        std::cout << "\nBest Energy: " << fitness_array[elites_index] << "\n";
        std::cout << "Best Curve: " << population.col(elites_index) << "\n";
    }
    

} // end of namespace
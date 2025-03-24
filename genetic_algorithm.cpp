#include <genetic_algorithm.h>
#include <vector>

namespace GA{

    // Constructor
    GeneticAlgorithm::GeneticAlgorithm(std::vector<double> starting_curve, int population_size, int num_generations, double mutation_rate)
    : starting_curve(starting_curve), population_size(population_size), num_generations(num_generations), mutation_rate(mutation_rate) {}


    
    

} // end of namespace
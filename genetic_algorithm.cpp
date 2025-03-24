#include <genetic_algorithm.h>
#include <microstrip_arrangement.h>
#include <vector>

namespace GA{

    // Constructor
    GeneticAlgorithm::GeneticAlgorithm(const MSA::MicrostripArrangement& arrangement,
        std::vector<double> starting_curve, 
        int population_size, 
        int num_generations, 
        double mutation_rate)
        : 
        arrangement(arrangement),
        starting_curve(starting_curve), 
        population_size(population_size), 
        num_generations(num_generations), 
        mutation_rate(mutation_rate) {}


    
    

} // end of namespace
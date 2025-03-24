#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>
namespace GA{
    class GeneticAlgorithm {
        private:
            // GA Properties
            std::vector<double> starting_curve;
            const int population_size;
            const int num_generations;
            const double mutation_rate;

            // Microstrip arrangement properties
            const MSA::MicrostripArrangement arrangement;

            //double fitnessFunction(const std::vector<double>& individual) const;
            
        public:
            // Constructor
            GeneticAlgorithm(const MSA::MicrostripArrangement& arrangement, 
                std::vector<double> starting_curve, 
                int population_size, 
                int num_generations, 
                double mutation_rate);

    };

} // end of namespace
#endif
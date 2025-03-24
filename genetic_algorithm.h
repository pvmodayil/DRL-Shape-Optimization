#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>

namespace GA{
    class GeneticAlgorithm {
        private:
            // GA Properties
            std::vector<double> starting_curveY;
            std::vector<double> starting_curveX;
            const int population_size;
            const int num_generations;
            const double mutation_rate;

            // Microstrip arrangement properties
            const MSA::MicrostripArrangement arrangement;

            double calculateFitness(const std::vector<double>& individual) const;
            
        public:
            // Constructor
            GeneticAlgorithm(const MSA::MicrostripArrangement& arrangement, 
                std::vector<double>& starting_curveY,
                std::vector<double>& starting_curveX, 
                int population_size, 
                int num_generations, 
                double mutation_rate);

    };

} // end of namespace
#endif
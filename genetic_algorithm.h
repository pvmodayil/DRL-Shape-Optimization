#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
namespace GA{

    class GeneticAlgorithm {
        private:
            std::vector<double> starting_curve;
            const int population_size;
            const int num_generations;
            const double mutation_rate;

        public:
            // Constructor
            GeneticAlgorithm(std::vector<double> starting_curve, int population_size, int num_generations, double mutation_rate);
    };

} // end of namespace
#endif
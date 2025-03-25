#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>

namespace GA{
    class GeneticAlgorithm {
        private:
            // GA Properties
            const Eigen::ArrayXd starting_curveY;
            const Eigen::ArrayXd starting_curveX;
            const int population_size;
            const int num_generations;
            const double mutation_rate;
            const MSA::MicrostripArrangement arrangement;
            
            // Functions
            Eigen::MatrixXd initializePopulation(const double& noise_scale) const;

            double calculateFitness(Eigen::ArrayXd& individual) const;
            
        public:
            // Constructor
            GeneticAlgorithm(MSA::MicrostripArrangement& arrangement, 
                Eigen::ArrayXd& starting_curveY,
                Eigen::ArrayXd& starting_curveX, 
                int population_size, 
                int num_generations, 
                double mutation_rate);
            
            // Main function to run the optimization
            void optimize(const double& noise_scale);

    };

} // end of namespace
#endif
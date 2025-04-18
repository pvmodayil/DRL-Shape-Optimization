#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>

namespace GA{
    class GeneticAlgorithm {
        private:
            // GA Properties
            // ------------------------------------------------------
            Eigen::ArrayXd starting_curveY;
            Eigen::ArrayXd starting_curveX;
            int population_size;
            int num_generations;
            double mutation_rate;
            MSA::MicrostripArrangement arrangement;
            
            // Functions
            // ------------------------------------------------------
            Eigen::MatrixXd initializePopulation(double& noise_scale);

            double calculateFitness(Eigen::ArrayXd& individual);
            
            // Parent selection  
            std::vector<size_t> selectElites(const Eigen::ArrayXd& fitness_array);
            std::vector<size_t> selectParents(const std::vector<size_t>& elites_indices, const Eigen::ArrayXd& fitness_array);

            // Reproduction
            Eigen::MatrixXd reproduce(Eigen::MatrixXd& population, Eigen::ArrayXd& fitness_array, double& noise_scale, std::vector<double>& energies);
            Eigen::VectorXd crossover(Eigen::VectorXd& parent1, Eigen::VectorXd& parent2);
            
        public:
            // Constructor
            GeneticAlgorithm(MSA::MicrostripArrangement& arrangement, 
                Eigen::ArrayXd& starting_curveY,
                Eigen::ArrayXd& starting_curveX, 
                int population_size, 
                int num_generations, 
                double mutation_rate);
            
            // Main function to run the optimization
            void optimize(double& noise_scale);

    };

} // end of namespace
#endif
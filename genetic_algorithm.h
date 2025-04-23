#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>
#include <tuple>

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
            size_t selectParent(const Eigen::ArrayXd& fitness_array);

            // Reproduction
            Eigen::MatrixXd reproduce(Eigen::MatrixXd& population, Eigen::ArrayXd& fitness_array, double& noise_scale);
            std::tuple<Eigen::VectorXd, Eigen::VectorXd> crossover(Eigen::VectorXd& parent1, Eigen::VectorXd& parent2, double eta=1.5);
            
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
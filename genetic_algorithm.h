#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <microstrip_arrangement.h>
#include <vector>
#include <tuple>
#include <random>
#include <omp.h>
namespace GA{
    // Result struct
    struct GeneticResult {
        Eigen::ArrayXd energy_convergence;
        Eigen::ArrayXd best_curve;
        double best_energy;

        GeneticResult(Eigen::ArrayXd energy_convergence, Eigen::VectorXd best_curve, double best_energy)
            : energy_convergence(energy_convergence), best_curve(best_curve), best_energy(best_energy) {}
    };
    
    // Genetic Algorithm Class
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

            // Random device
            std::vector<std::mt19937> rng_engines; // Vector of random number generators for parallel processing for thread safety
            //std::mt19937 rng;
            std::uniform_int_distribution<> parent_index_dist;
            
            // Functions
            // ------------------------------------------------------
            Eigen::MatrixXd initializePopulation(double& noise_scale);

            double calculateFitness(Eigen::ArrayXd& individual);
            
            // Parent selection  
            size_t selectElites(const Eigen::ArrayXd& fitness_array);
            size_t selectParent(const Eigen::ArrayXd& fitness_array, const int& thread_id);

            // Reproduction
            Eigen::MatrixXd reproduce(Eigen::MatrixXd& population, Eigen::ArrayXd& fitness_array, double& noise_scale);
            void crossover(Eigen::VectorXd& parent1, 
                Eigen::VectorXd& parent2, 
                Eigen::Ref<Eigen::VectorXd> child1, 
                Eigen::Ref<Eigen::VectorXd> child2, 
                double eta=1.5);
            
        public:
            // Constructor
            GeneticAlgorithm(MSA::MicrostripArrangement& arrangement, 
                Eigen::ArrayXd& starting_curveY,
                Eigen::ArrayXd& starting_curveX, 
                int population_size, 
                int num_generations, 
                double mutation_rate);
            
            // Main function to run the optimization
            void optimize(double& noise_scale, GeneticResult& result);

    };

} // end of namespace
#endif
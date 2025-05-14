#include "file_io.h"
#include "microstrip_arrangement.h"
#include "genetic_algorithm.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <omp.h>

// External Includes
#include <Eigen/Dense>
#include <pybind11/pybind11.h>

/*
*******************************************************
*            Shape Optimizer Function                 *
*******************************************************
*/
GA::GeneticResult optimizeShape(const double V0, 
        const double hw_micrstr, // Half-width of microstrip in meters 
        const double ht_micrstr, // Height of microstrip in meters 
        const double hw_arra, // Half-width of array in meters 
        const double ht_arra, // Height of array in meters 
        const double ht_subs, // Height of substrate in meters 
        const double er1, // Relative permittivity of air
        const double er2, // Relative permittivity of substrate
        const int N, // Number of Fourier coefficients
        const int num_pts, // Number of points in the bezier curve
        const Eigen::ArrayXd& action) { // Action vector
    
    // Create the microstrip arrangement
    //------------------------------------------------------
    MSA::MicrostripArrangement arrangement(V0, hw_micrstr, ht_micrstr, hw_arra, ht_arra, ht_subs, er1, er2, N);
    
    // Create the bezier curve
    //------------------------------------------------------
    Eigen::ArrayXd gptsX(num_pts);
    Eigen::ArrayXd gptsY(num_pts);
    MSA::getBezierCurve(action, hw_micrstr, hw_arra, num_pts, gptsX, gptsY);
    // Check if the curve is monotonically decreasing
    if (!MSA::isMonotonicallyDecreasing(gptsY)) {
        throw std::runtime_error("The curve is not monotonically decreasing.");
        return;
    }

    // Preprocess the curve points (the vn calculation takes the first and last points into account)
    //------------------------------------------------------
    Eigen::ArrayXd gptsX_filtered = gptsX.segment(1, num_pts - 2);
    Eigen::ArrayXd gptsY_filtered = gptsY.segment(1, num_pts - 2);

    // Initial energy calculation
    Eigen::ArrayXd vn = MSA::calculatePotentialCoeffs(arrangement.V0,
        arrangement.hw_micrstr,
        arrangement.hw_arra,
        arrangement.N,
        gptsY_filtered,
        gptsX_filtered);
    double init_energy = MSA::calculateEnergy(arrangement.er1,
        arrangement.er2,
        arrangement.hw_arra,
        arrangement.ht_arra,
        arrangement.ht_subs,
        arrangement.N,
        vn);
    
    // Genetic Algorithm class
    //------------------------------------------------------
    int population_size = 100; // Number of individuals in the population
    int num_generations = 1000; // Number of generations to evolve
    GA::GeneticAlgorithm ga_problem = GA::GeneticAlgorithm(arrangement,gptsY_filtered,gptsX_filtered,population_size,num_generations);
    double noise_scale = 0.1;       
    // Result struct
    double execution_time = 0.0;
    GA::GeneticResult result = GA::GeneticResult(Eigen::ArrayXd::Zero(num_generations+1), 
                                                Eigen::VectorXd(num_pts), 
                                                init_energy,
                                                execution_time); // Empty initialization
    result.energy_convergence(0) = init_energy; // Initial energy
    
    // Optimize
    //------------------------------------------------------
    auto start = std::chrono::high_resolution_clock::now();
    ga_problem.optimize(noise_scale, result);
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    std::chrono::duration<double, std::milli> duration = end - start;
    result.execution_time = duration.count(); // Convert to seconds

    return result;
}

int main(){
    return 0;
}
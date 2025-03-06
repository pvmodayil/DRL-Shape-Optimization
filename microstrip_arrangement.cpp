#include <iostream>
#include <format>
#include <string>
#include <utility>
#include <cassert>
#include <stdexcept>
#include <vector>

/*
*******************************************************
*            Necessary Conditons Chheck               *
*******************************************************
*/
bool isMonotonicallyDecreasing(const std::vector<double>& g){
    // Check if the vector has less than 2 elements
    if (g.size() < 2) {
        return true; // A single element or empty array is considered monotonically decreasing
    }

    // Iterate through the array and check if each element is greater than the next one
    for (size_t idx = 1; idx < g.size(); ++idx) {
        if (g[idx] >= g[idx - 1]) {
            return false; // Found an instance where the order is not decreasing
        }
    }

    return true; // All checks passed, so it is monotonically decreasing
}

bool isConvex(const std::vector<double>& g) {
    // Check if the vector has at least 3 elements
    if (g.size() < 3) {
        return false; // Not enough points to determine convexity
    }

    // Calculate the second differences
    std::vector<double> dx2(g.size() - 2);
    for (size_t idx = 1; idx < g.size() - 1; ++idx) {
        dx2[idx - 1] = g[idx + 1] - 2 * g[idx] + g[idx - 1];
    }

    // Check if all second differences are non-negative
    for (const auto& diff : dx2) {
        if (diff < 0) {
            return false; // Found a negative second difference, not convex
        }
    }

    return true; // All checks passed, so it is convex
}

/*
*******************************************************
*            Potential & Potential Coeffs             *
*******************************************************
*/
// Filter Vectors to be within dimensions
void filterVectors(double hw_micrstr, 
                    double hw_arra, 
                    std::vector<double> g, 
                    std::vector<double> x,
                    std::vector<double>& filtered_x,
                    std::vector<double>& filtered_g){
    std::cout<<"Input vectors are being filtered\n";
    
    // Assuming that both x and g have the same dimensions
    std::string error_message = std::format("Dimensions of x-axis vector and g-points vector do not match!");
    assert((std::format("Dimensions of x-axis vector and g-point vector do not match!g: {}, x: {}",g.size(),x.size()), g.size() == x.size()));

    // Create clear and add to the filtered vectors what is required
    filtered_g.clear();
    filtered_x.clear();
    for(size_t idx=0; idx<x.size(); ++idx){
        if(x[idx] > hw_micrstr && x[idx] < hw_arra){
            filtered_x.push_back(x[idx]);
            filtered_g.push_back(g[idx]);
        }
    }
}

//std::tie(filtered_g, filtered_x), std::pair<std::vector<double>, std::vector<double>> {filtered_g,filtered_x}

// Calculate the potential coefficients
std::vector<double> calculatePotentialCoeffs(double V0,
                                    double hw_micrstr, 
                                    double hw_arra, 
                                    int num_fs, 
                                    std::vector<double> g, 
                                    std::vector<double> x){
                                    
    if(g.size() != x.size()){
        throw std::runtime_error(std::format("Dimensions of x-axis vector and g-point vector do not match!g: {}, x: {}",g.size(),x.size()));
    }

    if(x[0] <= hw_micrstr || x[x.size() - 1] >= hw_arra){
        // first pair of g,x are passed as value and the second is passed as reference so the original vectors itself will be filtered
        filterVectors(hw_micrstr,hw_arra,g,x,g,x);
    }

    
}
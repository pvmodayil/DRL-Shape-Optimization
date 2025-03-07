#include <iostream>
#include <tuple>
#include <format>
#include <string>
#include <algorithm>
#include <numeric>
#include <numbers>
#include <utility>
#include <cassert>
#include <stdexcept>
#include <vector>

// External Includes
#include <Eigen/Dense>

// Constants
constexpr double PI = std::numbers::pi;

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
                    const std::vector<double> original_g, 
                    const std::vector<double> original_x,
                    std::vector<double>& filtered_x,
                    std::vector<double>& filtered_g){
    std::cout<<"Input vectors are being filtered\n";
    
    // Assuming that both x and g have the same dimensions
    std::string error_message = std::format("Dimensions of x-axis vector and g-points vector do not match!");
    assert((std::format("Dimensions of x-axis vector and g-point vector do not match!g: {}, x: {}",
            original_g.size(),original_x.size()), 
                original_g.size() == original_x.size()));

    // Create clear and add to the filtered vectors what is required
    filtered_g.clear();
    filtered_x.clear();
    for(size_t idx=0; idx<original_x.size(); ++idx){
        if(original_x[idx] > hw_micrstr && original_x[idx] < hw_arra){
            filtered_x.push_back(original_x[idx]);
            filtered_g.push_back(original_g[idx]);
        }
    }
}

//std::tie(filtered_g, filtered_x), std::pair<std::vector<double>, std::vector<double>> {filtered_g,filtered_x}

// Calculate the potential coefficients
Eigen::ArrayXd calculatePotentialCoeffs(const double V0,
                                    const double hw_micrstr, 
                                    const double hw_arra, 
                                    const int num_fs, 
                                    std::vector<double>& g, 
                                    std::vector<double>& x){

    // Check if the values passed are within the allowed range else filter
    if(x[0] <= hw_micrstr || x.back() >= hw_arra){
        // first pair of g,x are passed as value and the second is passed as reference so the original vectors itself will be filtered
        filterVectors(hw_micrstr,hw_arra,g,x,g,x);
    }
    
    // Set m as the size of the input vectors after filtering if required
    size_t m = x.size();
    
    // Assert the requirements
    assert(("Not enough spline knots, at least two points are required between half width of microstrip and half width of arrangement!", m > 1));
    assert((std::format("Dimensions of x-axis vector and g-point vector do not match!g: {}, x: {}",g.size(),x.size()), g.size() == x.size()));

    // For the next steps it is a must that the vectors are in the correct order
    if (!std::is_sorted(x.begin(), x.end())) {
        throw std::runtime_error("Input x-axis values are not sorted.");
    }

    // Calculate the potential coefficients
    // Create a array of Fourier coefficients
    Eigen::ArrayXd n = Eigen::ArrayXd::LinSpaced(num_fs, 0, num_fs - 1); // Create an array from 0 to num_fs-1

    
    Eigen::ArrayXd outer_coeff = (2.0/hw_arra)*V0*(1.0 / ((2 * n + 1) * PI / (2.0 * hw_arra)).square());
    
    // Calculate v_n1 and v_n3 for all n
    Eigen::ArrayXd v_n1 = (g[0] - V0) / (x[0] - hw_micrstr) *
        ((2 * n + 1).cast<double>() * x[0] * PI / (2 * hw_arra)).cos() -
        ((2 * n.array() + 1).cast<double>() * hw_micrstr * PI / (2 * hw_arra)).cos();

    Eigen::ArrayXd v_n3 = (g.back()) / (hw_arra - x.back()) *
        ((2 * n + 1).cast<double>() * x.back() * PI / (2 * hw_arra)).cos();

    // Reshape x into a column vector (m x 1)
    Eigen::MatrixXd x_t(m, 1);
    for (size_t i = 0; i < m; ++i) {
        x_t(i, 0) = x[i];
    }

    // Convert the x and g vectors to arrays
    Eigen::ArrayXd x_array = Eigen::Map<const Eigen::ArrayXd>(x.data(), x.size());
    Eigen::ArrayXd g_array = Eigen::Map<const Eigen::ArrayXd>(g.data(), g.size());

    // Calculate cos1 and cos2
    // Xi => last m-1 values
    Eigen::ArrayXd cos1 = (x_t.bottomRows(m - 1).array() * PI / (2 * hw_arra) * (2 * n + 1)).cos(); // m-1x1 * 1xN = m-1xN
    // Xi-1 => first m-1 values
    Eigen::ArrayXd cos2 = (x_t.topRows(m - 1).array() * PI / (2 * hw_arra) * (2 * n + 1)).cos(); // m-1x1 * 1xN = m-1xN

    // Calculate fac1: g[1:m] - g[0:m-1]     
    Eigen::ArrayXd fac1 = (g_array.segment(1, m-1) - g_array.segment(0, m-1)) /
                        (x_array.segment(1, m-1) - x_array.segment(0, m-1)); // 1xm-1

    // Calculate v_n2: fac1 multiplied by the difference of cosines
    Eigen::ArrayXd v_n2 = fac1 * (cos1 - cos2); // 1xm-1 * m-1xN = 1xN

    Eigen::ArrayXd vn = outer_coeff*(v_n1+v_n2+v_n3); // 1xN

    return vn;
}

std::vector<double> calculatePotential(double hw_arra, int num_fs, Eigen::ArrayXd vn, std::vector<double> x){
    
}
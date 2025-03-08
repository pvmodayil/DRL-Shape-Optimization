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
constexpr double e0 = 8.854E-12;
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
                                    const int N, 
                                    std::vector<double>& g, 
                                    std::vector<double>& x){

    // Check if the values passed are within the allowed range else filter
    if(x[0] <= hw_micrstr || x.back() >= hw_arra){
        // first pair of g,x are passed as value and the second is passed as reference so the original vectors itself will be filtered
        filterVectors(hw_micrstr,hw_arra,g,x,g,x);
    }
    
    // Set M as the size of the input vectors after filtering if required
    size_t M = x.size();
    
    // Assert the requirements
    assert(M > 1 && "Not enough spline knots, at least two points are required between half width of microstrip and half width of arrangement!");
    if (g.size() != x.size()) {
        throw std::invalid_argument(std::format("Dimensions of x-axis vector and g-point vector do not match! g: {}, x: {}", g.size(), x.size()));
    }

    // For the next steps it is a must that the vectors are in the correct order
    if (!std::is_sorted(x.begin(), x.end())) {
        throw std::runtime_error("Input x-axis values are not sorted.");
    }

    // Create a array of Fourier coefficients
    Eigen::ArrayXd n = (Eigen::ArrayXd::LinSpaced(N, 0, N - 1)); // 1xN

    // Calculate outer coefficient
    Eigen::ArrayXd outer_coeff = (2.0/hw_arra)*V0*(1.0 / ((2 * n + 1) * PI / (2.0 * hw_arra)).square()); // 1xN
    
    // Calculate v_n1 and v_n3 for all n
    Eigen::ArrayXd v_n1 = (g[0] - V0) / (x[0] - hw_micrstr) *
        ((2 * n + 1) * x[0] * PI / (2 * hw_arra)).cos() -
        ((2 * n + 1) * hw_micrstr * PI / (2 * hw_arra)).cos(); // 1xN

    Eigen::ArrayXd v_n3 = (g[M]) / (hw_arra - x[M]) *
        ((2 * n + 1) * x[M] * PI / (2 * hw_arra)).cos(); // 1xN

    // Edge case when the input vector is very small
    if(M == 2){
        // Calculate cos1 and cos2 (segment takes start index and number of positions including start index to be taken)
        Eigen::ArrayXd cos1 = ((x[1] * PI / (2 * hw_arra)) * (2 * n + 1)).cos(); // 1xN 
        // Require values from first element to the second last element (0 to (M-1th))
        Eigen::ArrayXd cos2 = ((x[0] * PI / (2 * hw_arra)) * (2 * n + 1)).cos(); // 1xN

        Eigen::ArrayXd v_n2 = (g[1]-g[0])/(x[1]-x[0])*(cos1-cos2); // 1xN

        Eigen::ArrayXd vn = outer_coeff*(v_n1+v_n2+v_n3); // 1xN

        return vn;
    }

    // Convert the x and g vectors to arrays
    Eigen::ArrayXd x_array = Eigen::Map<const Eigen::ArrayXd>(x.data(), x.size(), 1); // Mx1
    Eigen::ArrayXd g_array = Eigen::Map<const Eigen::ArrayXd>(g.data(), g.size(), 1); // Mx1

    // Calculate cos1 and cos2 (segment takes start index and number of positions including start index to be taken)
    // Require values from second element to the last element (1 to (M-1)th)
    Eigen::ArrayXd cos1 = ((x_array.bottomRows(M - 1) * PI / (2 * hw_arra)).matrix() * (2 * n + 1).matrix()).array().cos(); // M-1x1 * 1xN = M-1xN
    // Require values from first element to the second last element (0 to (M-1th))
    Eigen::ArrayXd cos2 = ((x_array.topRows(M - 1) * PI / (2 * hw_arra)).matrix() * (2 * n + 1).matrix()).array().cos(); // M-1x1 * 1xN = M-1xN

    // Calculate fac1: g[1:M] - g[0:M-1]     
    Eigen::ArrayXd coeff_vn2 = (g_array.bottomRows(M-1) - g_array.topRows(M-1)) /
                        (x_array.bottomRows(M-1) - x_array.topRows(M-1)); // M-1x1

    // Calculate v_n2: fac1 multiplied by the difference of cosines
    Eigen::ArrayXd v_n2 = coeff_vn2.transpose().matrix() * (cos1 - cos2).matrix(); // 1xM-1 * M-1xN = 1xN

    Eigen::ArrayXd vn = outer_coeff*(v_n1+v_n2+v_n3); // 1xN

    assert(vn.rows() == 1 && "The calculated potential coefficients resulted in a column vector");

    return vn;
}

Eigen::ArrayXd calculatePotential(const double hw_arra, 
                                    const int num_fs, 
                                    Eigen::ArrayXd& vn, 
                                    std::vector<double>& x){

    assert(vn.rows() != 0 && "Potenntial coefficients vn is empty"); // && message is a trick to print message in assert as the assert checks failure of conditons                                    
    
    // Create the Fourier coefficients
    Eigen::ArrayXd n = Eigen::ArrayXd::LinSpaced(num_fs, 0, num_fs - 1).transpose(); // Nx1

    // Convert the x vector to a MatrixXd (Mxm)
    Eigen::MatrixXd x_matrix = Eigen::Map<const Eigen::MatrixXd>(x.data(), 1, x.size()); // 1xM

    // Calculate cosines: NxM => the final .array() will make it an array type
    Eigen::ArrayXd cos1 = (((2 * n + 1) * (PI / (2 * hw_arra))).matrix() * x_matrix).array().cos(); // NxM

    // Multiply vn with cos1; vn should be a column vector for proper broadcasting
    // Expected dimension of vn is 1xN but verify it
    if(vn.rows() > 1){
        // vn is a column vector, convert it to row vector and calculate
        Eigen::ArrayXd VF = (vn.transpose().matrix() * cos1.matrix()).array(); // 1xN * NxM = 1xM
        return VF;
    }

    Eigen::ArrayXd VF = (vn.matrix() * cos1.matrix()).array(); // 1xN * NxM = 1xM

    return VF;
}

/*
*******************************************************
*                      Energy                         *
*******************************************************
*/

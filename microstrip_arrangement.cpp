#include "microstrip_arrangement.h"
#include <iostream>
#include <format>
#include <string>
#include <numbers>
#include <cassert>
#include <stdexcept>
#include <vector>

// External Includes
#include <Eigen/Dense>

// Constants
constexpr double PI = std::numbers::pi;
constexpr double c0 = 2.99792458e8; // Speed of light in vacuum in m/s

namespace MSA{

    /*
    *******************************************************
    *            Bezier Curve Generator                   *
    *******************************************************
    */
    void getBezierCurve(const Eigen::ArrayXd& action, 
        const double& hw_micrstr, 
        const double& hw_arra, 
        const int& num_pts,
        Eigen::ArrayXd& gptsX,
        Eigen::ArrayXd& gptsY) {
        
        // Control points
        Eigen::Vector2d P0(hw_micrstr, 1.0);
        double p10 = hw_micrstr + (hw_arra - hw_micrstr) * action(0);
        Eigen::Vector2d P1(p10, action(1));
        double p20 = hw_micrstr + (hw_arra - hw_micrstr) * action(2);
        Eigen::Vector2d P2(p20, action(3)); 
        Eigen::Vector2d P3(hw_arra, 0.0);
        // Stack the control points to do matrix multiplication
        Eigen::MatrixXd control_points(4, 2);
        control_points.row(0) = P0;
        control_points.row(1) = P1;
        control_points.row(2) = P2;
        control_points.row(3) = P3;
        
        // Generate the bezier curve coefficients
        Eigen::ArrayXd t = Eigen::ArrayXd::LinSpaced(num_pts, 0.0, 1.0); // 1xnum_pts
        Eigen::VectorXd B0 = (1 - t).pow(3);
        Eigen::VectorXd B1 = 3 * t * (1 - t).pow(2);
        Eigen::VectorXd B2 = 3 * t.pow(2) * (1 - t);
        Eigen::VectorXd B3 = t.pow(3);
        
        // Stack the coefficients to do matrix multiplication
        Eigen::MatrixXd curve_coefficients(num_pts, 4);
        curve_coefficients.col(0) = B0;
        curve_coefficients.col(1) = B1;
        curve_coefficients.col(2) = B2;
        curve_coefficients.col(3) = B3;

        // Calculate the curve points using matrix multiplication
        Eigen::MatrixXd curve_points = curve_coefficients * control_points;

        // Extract x and y coordinates
        gptsX = curve_points.col(0).array();
        gptsY = curve_points.col(1).array();
    }
    /*
    *******************************************************
    *            Necessary Conditons Check               *
    *******************************************************
    */
    bool isMonotonicallyDecreasing(Eigen::ArrayXd& g){
        size_t n = g.size();

        // Check if the vector has less than 2 elements
        if (n < 2) {
            return true; // A single element or empty array is considered monotonically decreasing
        }
        
        return (g.tail(n) <= g.head(n)).all();
    }

    bool isConvex(Eigen::ArrayXd& g){
        size_t n = g.size();

        // Check if the vector has at least 3 elemen
        if (n < 3) {
            return false; // Not enough points to determine convexity
        }

        // Calculate the second differences
        Eigen::ArrayXd dx2(n - 2);
        // segment(start idx, length)
        dx2 = g.segment(2, n-2) 
            - 2.0 * g.segment(1, n-2) 
            + g.segment(0, n-2);

        // Check if all second differences are non-negative
        return (dx2 >= 0).all(); 
    }

    /*
    *******************************************************
    *            Potential & Potential Coeffs             *
    *******************************************************
    */
    // Filter Vectors to be within dimensions
    void filterVectors(const double& hw_micrstr, 
                        const double& hw_arra, 
                        const std::vector<double> original_g, 
                        const std::vector<double> original_x,
                        std::vector<double>& filtered_g,
                        std::vector<double>& filtered_x){
        // Assuming that both x and g have the same dimensions
        std::cout<<"Input vectors are being filtered\n";

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

    //std::tie(filtered_g, filtered_x), std::pair<std::vector<double>, std::vector<double>> {filtered_g,filtered_x} #include <tuple>

    // Calculate the potential coefficients
    Eigen::ArrayXd calculatePotentialCoeffs(const double& V0,
                                        const double& hw_micrstr, 
                                        const double& hw_arra, 
                                        const int& N, 
                                        const Eigen::ArrayXd& g, 
                                        const Eigen::ArrayXd& x){

        // Set M as the size of the input vectors after filtering if required
        size_t M = x.size();
        
        // Assert the requirements
        if(M < 2){
            throw std::invalid_argument("Not enough spline knots, at least two points are required between half width of microstrip and half width of arrangement!");
        }
        
        if (g.size() != x.size()) {
            throw std::invalid_argument(std::format("Dimensions of x-axis vector and g-point vector do not match! g: {}, x: {}", g.size(), x.size()));
        }

        // For the next steps it is a must that the vectors are in the correct order
        if (!std::is_sorted(x.begin(), x.end())) {
            throw std::runtime_error("Input x-axis values are not sorted.");
        }

        // Create an array of Fourier coefficients
        Eigen::ArrayXd n = (Eigen::ArrayXd::LinSpaced(N, 0, N - 1)); // Nx1

        // Calculate outer coefficient
        Eigen::ArrayXd outer_coeff = (2.0/hw_arra) * V0 * (1.0/((2*n+1)*PI/(2.0*hw_arra)).square()); // Nx1

        // Calculate v_n1 and v_n3 for all n
        Eigen::ArrayXd v_n1 = ((g[0]-V0) / (x[0]-hw_micrstr)) * (
            ((2*n+1)*x[0]*PI/(2*hw_arra)).cos()-((2*n+1)*hw_micrstr*PI/(2*hw_arra)).cos()
        ); // Nx1
        
        Eigen::ArrayXd v_n3 = (g[M-1] / (hw_arra - x[M-1])) * (
            ((2*n+1) * x[M-1] * PI / (2*hw_arra)).cos()
        ); // Nx1
        
        // Edge case when the input vector is very small
        if(M == 2){
            // Calculate cos1 and cos2 (segment takes start index and number of positions including start index to be taken)
            Eigen::ArrayXd cos1 = ((x[1] * PI / (2 * hw_arra)) * (2 * n + 1)).cos(); // Nx1 
            // Require values from first element to the second last element (0 to (M-1th))
            Eigen::ArrayXd cos2 = ((x[0] * PI / (2 * hw_arra)) * (2 * n + 1)).cos(); // Nx1

            Eigen::ArrayXd v_n2 = (g[1]-g[0])/(x[1]-x[0])*(cos1-cos2); // Nx1

            Eigen::ArrayXd vn = outer_coeff*(v_n1+v_n2+v_n3); // Nx1

            return vn;
        }
        
        // Calculate cos1 and cos2 (segment takes start index and number of positions including start index to be taken)
        // Require values from second element to the last element (1 to (M-1)th)
        Eigen::ArrayXXd cos1 = ((x.bottomRows(M - 1) * PI / (2 * hw_arra)).matrix() * (2 * n + 1).matrix().transpose()).array().cos(); // M-1x1 * 1xN = M-1xN
        // Require values from first element to the second last element (0 to (M-1th))
        Eigen::ArrayXXd cos2 = ((x.topRows(M - 1) * PI / (2 * hw_arra)).matrix() * (2 * n + 1).matrix().transpose()).array().cos(); // M-1x1 * 1xN = M-1xN
        // Calculate fac1: g[1:M] - g[0:M-1]     
        Eigen::ArrayXd coeff_vn2 = (g.bottomRows(M-1) - g.topRows(M-1)) /
                            (x.bottomRows(M-1) - x.topRows(M-1)); // M-1x1

        // Calculate v_n2: fac1 multiplied by the difference of cosines
        Eigen::MatrixXd cos_diff = (cos1 - cos2).matrix(); // M-1xN
        Eigen::MatrixXd coeff_matrix = coeff_vn2.matrix(); // M-1x1
        Eigen::MatrixXd v_n2 = (coeff_matrix.transpose() * cos_diff); // 1xM-1 * M-1xN = 1xN
        
        Eigen::ArrayXd vn = outer_coeff*(v_n1+(v_n2.transpose().array())+v_n3); // Nx1

        return vn;
    }

    Eigen::ArrayXd calculatePotential(const double& hw_arra, 
                                    const int& N, 
                                    Eigen::ArrayXd& vn, // Expected to be Nx1 dimension
                                    std::vector<double>& x){

        // Input potential coefficients must not be empty
        if(vn.rows() == 0){
            throw std::invalid_argument("Potenntial coefficients vn is empty");
        }                                  
        
        // Create the Fourier coefficients
        Eigen::ArrayXd n = Eigen::ArrayXd::LinSpaced(N, 0, N - 1); // Nx1

        // Convert the x vector to a MatrixXd
        Eigen::MatrixXd x_matrix = Eigen::Map<const Eigen::MatrixXd>(x.data(), 1, x.size()); // 1xM

        // Calculate cosines: NxM => the final .array() will make it an array type
        Eigen::ArrayXXd cos1 = (((2 * n + 1) * (PI / (2 * hw_arra))).matrix() * x_matrix).array().cos(); // NxM
        
        // Multiply vn with cos1; input vn is expected to be Nx1 dimension
        // vn is a column vector, convert it to row vector and calculate
        Eigen::ArrayXd VF = (vn.matrix().transpose() * cos1.matrix()).transpose().array(); // 1xN * NxM = 1xM => will be Mx1 after transpose
        // Since expected result is a 1D array ArrayXd is enough and not ArrayXXd
        return VF;   
    }

    /*
    *******************************************************
    *                      Energy                         *
    *******************************************************
    */
    // Function to compute the natural logarithm of the sinh function
    Eigen::ArrayXd logsinh(const Eigen::ArrayXd& vector){
        // sinh = (e^x - e^(-x))/2 => ln(sinh(x)) = ln(e^x - e^(-x)) - ln(2)

        // Get absolute values
        Eigen::ArrayXd absolute_vector = vector.abs();

        // Exponential of the values after normalization to avoid overflow
        Eigen::ArrayXd exp_x = vector.exp();
        Eigen::ArrayXd exp_neg_x = (-vector).exp();
        Eigen::ArrayXd small_result = (exp_x - exp_neg_x).log() - std::log(2.0);
        Eigen::ArrayXd large_result = absolute_vector - std::log(2.0);
        
        // Select values based on conditions, large x values are scaled
        return (absolute_vector > 33.0).select(large_result, small_result);
    }

    // Function to compute the natural logarithm of the cosh function
    Eigen::ArrayXd logcosh(const Eigen::ArrayXd& vector){
        // cosh = (e^x + e^(-x))/2 => ln(cosh(x)) = ln(e^x + e^(-x)) - ln(2)

        // Get absolute values
        Eigen::ArrayXd absolute_vector = vector.abs();

        // Exponential of the values after normalization to avoid overflow
        Eigen::ArrayXd exp_x = vector.exp();
        Eigen::ArrayXd exp_neg_x = (-vector).exp();
        Eigen::ArrayXd small_result = (exp_x + exp_neg_x).log() - std::log(2.0);
        Eigen::ArrayXd large_result = absolute_vector - std::log(2.0);
        
        // Select values based on conditions, large x values are scaled
        return (absolute_vector > 33.0).select(large_result, small_result);
    }

    // Function to calculate the energy of the curve
    double calculateEnergy(const double& er1,
                        const double& er2,
                        const double& hw_arra,
                        const double& ht_arra,
                        const double& ht_subs,
                        const int& N,
                        Eigen::ArrayXd& vn){
        
        // Set the constants
        double e0 = 8.854E-12; // Permittivity constant e0
        double e1 = er1 * e0;
        double e2 = er2 * e0;

        // Create a array of Fourier coefficients
        Eigen::ArrayXd n = (Eigen::ArrayXd::LinSpaced(N, 0, N - 1)); // 1xN

        // Coefficients
        Eigen::ArrayXd coeff = (2 * n + 1) * PI * vn.square(); // 1xN

        // w1
        Eigen::ArrayXd logcosh1 = logcosh((2 * n + 1) * PI * (ht_arra - ht_subs) / (2 * hw_arra)); // 1xN
        Eigen::ArrayXd logsinh1 = logsinh((2 * n + 1) * PI * (ht_arra - ht_subs) / (2 * hw_arra)); // 1xN
        Eigen::ArrayXd coth1 = (logcosh1 - logsinh1).exp(); // 1xN
        double w1 = e1 / 4 * coeff.matrix().dot(coth1.matrix());

        // w2
        Eigen::ArrayXd logcosh2 = logcosh((2 * n + 1) * PI * ht_subs / (2 * hw_arra));
        Eigen::ArrayXd logsinh2 = logsinh((2 * n + 1) * PI * ht_subs / (2 * hw_arra));
        Eigen::ArrayXd coth2 = (logcosh2 - logsinh2).exp();
        double w2 = e2 / 4 * coeff.matrix().dot(coth2.matrix());

        double W12 = w1 + w2;

        return W12;
    }

    void calculateImpedance(const double& V0, 
        const double& WD, 
        const double& WL, 
        double& ZD, 
        double& ZL,
        double& epsilon_eff){
        // Calculate the impedance of the microstrip line
        double CD = 2.0 * WD/(V0*V0); // Capacitance of the microstrip line case D
        double CL = 2.0 * WL/(V0*V0); // Capacitance of the microstrip line case L

        ZD = 1.0 / (c0 * std::sqrt(CD*CL)); // Impedance of the microstrip line case D
        ZL = 1.0 / (c0 * CL); // Impedance of the microstrip line case L

        epsilon_eff = CD / CL; // Effective permittivity of the microstrip line
        }
} // namespace MSA
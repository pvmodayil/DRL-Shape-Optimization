#ifndef MICROSTRIP_ARRANGEMET_H
#define MICROSTRIP_ARRANGEMET_H

#include <format>
#include <string>
#include <numbers>
#include <cassert>
#include <stdexcept>
#include <vector>

// External Includes
#include <Eigen/Dense>

namespace MSA{

bool isMonotonicallyDecreasing(const std::vector<double>& g);

bool isConvex(const std::vector<double>& g);

void filterVectors(const double& hw_micrstr, 
    const double& hw_arra, 
    const std::vector<double> original_g, 
    const std::vector<double> original_x,
    std::vector<double>& filtered_x,
    std::vector<double>& filtered_g);
Eigen::ArrayXd calculatePotentialCoeffs(const double& V0,
    const double& hw_micrstr, 
    const double& hw_arra, 
    const int& N, 
    std::vector<double> g, 
    std::vector<double> x);

Eigen::ArrayXd calculatePotential(const double& hw_arra, 
    const int& N, 
    Eigen::ArrayXd& vn, 
    std::vector<double>& x);

Eigen::ArrayXd logsinh(const Eigen::ArrayXd& vector);

Eigen::ArrayXd logcosh(const Eigen::ArrayXd& vector);

double calculateEnergy(const double& er1,
    const double& er2,
    const double& hw_arra,
    const double& ht_arra,
    const double& ht_subs,
    const int& N,
    Eigen::ArrayXd& vn);

} // namespace MSA

#endif
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

struct MicrostripArrangement {
    const double V0;
    const double hw_micrstr;
    const double ht_micrstr;
    const double hw_arra;
    const double ht_arra;
    const double ht_subs;
    const double er1;
    const double er2;
    const int N; 

    MicrostripArrangement(
        const double V0, 
        const double hw_micrstr, // Half-width of microstrip in meters 
        const double ht_micrstr, // Height of microstrip in meters 
        const double hw_arra, // Half-width of array in meters 
        const double ht_arra, // Height of array in meters 
        const double ht_subs, // Height of substrate in meters 
        const double er1, // Relative permittivity of air
        const double er2, // Relative permittivity of substrate
        const int N) // Number of Fourier coefficients
        : 
        V0(V0), 
        hw_micrstr(hw_micrstr), 
        ht_micrstr(ht_micrstr), 
        hw_arra(hw_arra), 
        ht_arra(ht_arra), 
        ht_subs(ht_subs), 
        er1(er1), 
        er2(er2), 
        N(N){}
};

bool isMonotonicallyDecreasing(Eigen::ArrayXd& g);

bool isConvex(Eigen::ArrayXd& g);

void filterVectors(const double& hw_micrstr, 
    const double& hw_arra, 
    const std::vector<double> original_g, 
    const std::vector<double> original_x,
    std::vector<double>& filtered_g,
    std::vector<double>& filtered_x);
Eigen::ArrayXd calculatePotentialCoeffs(const double& V0,
    const double& hw_micrstr, 
    const double& hw_arra, 
    const int& N, 
    const Eigen::ArrayXd& g, 
    const Eigen::ArrayXd& x);

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
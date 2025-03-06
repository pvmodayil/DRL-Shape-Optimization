#include <iostream>
#include <vector>

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

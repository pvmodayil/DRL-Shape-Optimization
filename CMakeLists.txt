cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ShapeOptimization)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable compiler optimizations (e.g., -O3 for maximum optimization)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Specify the executable target
add_executable(microstrip_arrangement main.cpp file_io.cpp microstrip_arrangement.cpp genetic_algorithm.cpp)

# Include directories
include_directories(.)

# Enable OpenMP if available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(microstrip_arrangement PUBLIC OpenMP::OpenMP_CXX)
endif()
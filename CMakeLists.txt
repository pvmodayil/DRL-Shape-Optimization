cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ShapeOptimization)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the executable target
add_executable(microstrip_arrangement main.cpp file_io.cpp microstrip_arrangement.cpp)

# Include directories
include_directories(.)
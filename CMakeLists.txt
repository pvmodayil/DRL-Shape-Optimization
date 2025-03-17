cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ShapeOptimization)

# Specify the executable target
add_executable(microstrip_arrangement main.cpp microstrip_arrangement.cpp)

# Include directories
include_directories(.)
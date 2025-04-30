# DRL-Shape-Optimization

## Overview
This project implements a deep reinforcement learning (DRL) agent, specifically using the Soft Actor-Critic (SAC) algorithm, to predict the shape of a potential curve across a PCB single microstrip line based on various PCB parameters. The initial predicted shape is further optimised using a Genetic Algorithm (GA).

## Background
The potential curve shape is determined using Thomson's theorem, which states that the potential curve of the microstrip arrangement will be the one with the least energy. This principle serves as the foundation for our optimisation approach.

## Features
- **Deep Reinforcement Learning**: Utilises SAC to predict initial shapes.
- **Genetic Algorithm Optimization**: Implements GA in C++ for enhanced control and faster execution speeds.
  
## Getting Started
To get started with this project, clone the repository and follow these steps:

1. Clone this repository:
   ```bash
   git clone https://github.com/pvmodayil/DRL-Shape-Optimization.git
   cd DRL-Shape-Optimization
   ```
2. Build and run
    ```bash
    mkdir build
    cd build
    cmake ..
    cmake --build .
    .\microstrip_arrangement.exe
    ```
## Further information
For more detailed information about this study, please take a look at my thesis: [ShapeOptimization-Thesis](https://github.com/pvmodayil/MasterThesis-Shape-Optimisation-Using-DRL).


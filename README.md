[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/tKSbaXxd)
# challenge3 - 2D Laplace Equation Solver (Hybrid MPI + OpenMP)

## Goal of the Challenge
The goal of this project is to design and implement a high-performance solver for the 2D Laplace equation. 
The code uses the Finite Difference Method (FDM) and adopts a **hybrid** parallelization strategy:
* **MPI (Message Passing Interface):** Used for 1D row-wise domain decomposition and boundary exchange (ghost cells) between processes in distributed memory.
* **OpenMP:** Used to accelerate the internal computations within individual domains in shared memory by leveraging multi-threading.

Moreover, the project includes an advanced implementation of the **Schwarz Alternating Method** (Block Jacobi with overlap) to improve convergence speed, togheter with tools for automated performance analysis, plot generation, and 3D output for ParaView.

## Project Structure
The repository is organized to separate interfaces, source files and the testing phase:
* `Makefile`: automates the project build process
* `README.md`: this documentation file
* `Doxyfile`: configuration file for Doxygen documentation generation
* `include`: folder with all the header files
    - `Functions.hpp`: Defines analytical functions and test cases (forcing term, BCs)
    - `Matrix.hpp`: Custom class to handle 2D matrices in 1D arrays (optimized for MPI)
    - `MpiDomain.hpp`: Struct for automatic load-balancing and domain decomposition
    - `ParallelSchwarzSolver.hpp`: Interface for Schwarz method
    - `ParallelSolver.hpp`: Interface for the standard parallel solver (MPI + OpenMP)
    - `SerialSolver.hpp`: Interface for the serial solver
    - `VTKExport.hpp`: Interface for exporting results to 3D format
* `src`: folder with all the source file, with the implementations of the aforementioned .hpp files
    - `main.cpp`: Parses arguments, runs solvers and record timings for a first test comparing serial and parallel Jacobi solver
    - `Functions.cpp`: Implementation of exact test cases
    - `ParallelSchwarzSolver.cpp`: Logic for the Schwarz solver
    - `ParallelSolver.cpp`: Logic for the parallel Jacobi solver
    - `SerialSolver.cpp`: Logic for the serial solver
    - `VTKExport.cpp`: Logic to gather domains and write .vtk files
* `test`: automated testing environment
    - `data`: Auto-generated folder to store CSV results
    - `output`: Auto-generated folder that stores 3 plots: execution times, speedup and convergence error
    - `plot_results.py`
    - `run_tests.sh`: Bash script to run scalability tests across grids and processes
* `docs`: (Generated) HTML documentation created by Doxygen

# Compilation and Execution
# Prerequisites
- C++ Compiler (C++17 support)
- OpenMPI (or equivalent MPI implementation) and OpenMP
- Python 3 (with pandas and matplotlib libraries for plotting)
# Compilation
To compile the code, navigate to the root directory and run: `make`
Then, to run the program manually by specifying the grid size (e.g., N=256) and the number of processes (e.g., 4): `mpirun -np 4 ./main.out 256`
# Automated testing and plot generation
The project includes an automated performance analysis tool. To generate the plots please do the following commands:
* Running the test script: The following script executes the program across different grid sizes using 1, 2, and 4 MPI processes, saving execution times and errors into test/data/scalability_results.csv.
* Generating the Plots: By running the Python script, three .png images will be automatically generated in the grafici/ folder: `python3 plot_results.py`
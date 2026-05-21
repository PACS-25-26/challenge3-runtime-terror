#define OMPI_SKIP_MPICXX 1 // Avoid conflicts with MPI C++ bindings
#include "SerialSolver.hpp"
#include "ParallelSolver.hpp"
#include "Functions.hpp"
#include <iostream>
#include <string>
#include <mpi.h>

int main(int argc, char* argv[]) {
    // Initial checks for command line arguments
    MPI_Init(&argc, &argv); // Initialize MPI environment

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processes

    if (argc != 2) {
        if (rank == 0) { // Only the master process should print the error message
            std::cerr << "Error: incorrect number of arguments." << std::endl;
        }
        MPI_Finalize(); // Clean up MPI environment
        return 1;
    }
    
    int n = std::stoi(argv[1]);
    
    if (n <= 2) {
        std::cerr << "Error: the dimension n must be > 2." << std::endl;
        MPI_Finalize(); // Clean up MPI environment
        return 1;
    }
    
    MpiDomain mpi_domain(n, rank, size); // Create the MPI domain configuration
    ParallelSolver solver(mpi_domain, EquationData::forcing_term);
    
    if (rank == 0) { // Only the master process should print this message
        std::cout << "Initialization completed on master process." << std::endl;
    }
    
    solver.solve(50000, 1e-6);
    
    MPI_Finalize(); // Clean up MPI environment
    return 0;
}
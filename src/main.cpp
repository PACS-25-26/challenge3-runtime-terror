#include "SerialSolver.hpp"
#include "ParallelSolver.hpp"
#include "ParallelSchwarzSolver.hpp"
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
    
    int max_iter = 50000;
    double tol = 1e-6;

    double time_serial = 0.0; // To store the serial execution time for speedup calculation
    double err_serial = 0.0; // To store the serial error for comparison
    double time_parallel = 0.0; // To store the parallel execution time for speedup calculation
    double err_parallel = 0.0; // To store the parallel error for comparison
    double time_schwarz = 0.0; // To store the Schwarz execution time for speedup calculation
    double err_schwarz = 0.0; // To store the Schwarz error for comparison

    // set up problem data (forcing term, exact solution, boundary condition)
    ProblemData problem = TestCases::Sine(); 

    // SERIAL SOLVER TEST (only on rank 0)
    if (rank == 0) {
        std::cout << "\nTEST GRIGLIA " << n << "x" << n << std::endl;
        std::cout << "Starting serial test..." << std::endl;
        
        SerialSolver serial_solver(n, problem.f, problem.g);
        
        double start_serial = MPI_Wtime();
        serial_solver.solve(max_iter, tol);
        double end_serial = MPI_Wtime();

        time_serial = end_serial - start_serial;
        
        err_serial = serial_solver.compute_analytical_error(problem.exact_solution);
        
        std::cout << "[SERIAL] Time: " << time_serial 
                  << " s | Error: " << err_serial << std::endl;
    }

    // PARALLEL SOLVER TEST - JACOBI 
    MPI_Barrier(MPI_COMM_WORLD);  // To start with a clean slate, we synchronize all processes here

    if (rank == 0) {
        std::cout << "\nStarting parallel test with " << size << " processes..." << std::endl;
    }

    MpiDomain domain(n, rank, size);
    ParallelSolver parallel_solver(domain, problem.f, problem.g);
    
    double start_parallel = MPI_Wtime();
    parallel_solver.solve(max_iter, tol);
    double end_parallel = MPI_Wtime();

    time_parallel = end_parallel - start_parallel;
    
    err_parallel = parallel_solver.compute_analytical_error(problem.exact_solution);

    // Only the master process prints the parallel results to avoid cluttering the output
    if (rank == 0) {
        std::cout << "[PARALLEL] Time: " << time_parallel 
                  << " s | Error: " << err_parallel << std::endl;
        std::cout << "\n Speedup: " << time_serial / time_parallel << std::endl;
    }
    
    parallel_solver.export_vtk("test/output/vtk/solution_jacobi.vtk");


    // PARALLEL SOLVER TEST - BLOCK JACOBI
    MPI_Barrier(MPI_COMM_WORLD); // To start with a clean slate, we synchronize all processes here
    if (rank == 0){
        std::cout << "\nStarting Schwarz test with " << size << " processes..." << std::endl;
    }

    int max_global_iter = 2000;
    int max_local_iter = 100;
    ParallelSchwarzSolver schwarz_solver(domain, problem.f, problem.g);

    double start_schwarz = MPI_Wtime();
    schwarz_solver.solve(max_global_iter, max_local_iter, tol);
    double end_schwarz = MPI_Wtime();

    time_schwarz = end_schwarz - start_schwarz;

    err_schwarz = schwarz_solver.compute_analytical_error(problem.exact_solution);

    if(rank == 0){
        std::cout << "[SCHWARZ] Time: " << time_schwarz 
                  << " s | Error: " << err_schwarz << std::endl;
        std::cout << "\n Speedup: " << time_serial / time_schwarz << std::endl;
    }

    schwarz_solver.export_vtk("test/output/vtk/solution_schwarz.vtk");

    MPI_Finalize(); // Clean up MPI environment
    return 0;
}
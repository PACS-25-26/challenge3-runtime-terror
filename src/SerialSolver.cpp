#include "SerialSolver.hpp"
#include <cmath>
#include <iostream>

// Constructor
SerialSolver::SerialSolver(int grid_size, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term)
    : n(grid_size), 
      h(1.0 / (grid_size - 1.0)), 
      U(grid_size, grid_size),       
      U_new(grid_size, grid_size),   
      f(forcing_term),
      bc(boundary_term) 
{}

void SerialSolver::solve(int max_iter, double tol) {
    double error = tol + 1.0; 
    int iter = 0;
    
    std::cout << "Initializing Serial Jacobi..." << std::endl;

    while (error > tol && iter < max_iter) {
        error = 0.0;

        // Only on internal nodes, the borders are fixed at 0.0
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                // Physical coordinates of the node (i, j)
                double x = i * h;
                double y = j * h;

                U_new(i, j) = 0.25 * (U(i - 1, j) + U(i + 1, j) + U(i, j - 1) + U(i, j + 1) + (h * h * f(x, y)));

                double diff = U_new(i, j) - U(i, j);
                error += diff * diff;
            }
        }

        // BCs
        for (int i = 0; i < n; ++i){
            for (int j = 0; j < n; ++j){
                if(i == 0 || i == n-1 || j == 0 || j == n-1){
                    double x = i * h;
                    double y = j * h;
                    U(i, j) = bc(x, y);
                    U_new(i, j) = bc(x, y);
                }
            }
        }

        error = std::sqrt(h * error);

        U.swap(U_new);
        
        iter++;
        
        // Print progress every 1000 iterations
        if (iter % 1000 == 0) {
            std::cout << "Iteration " << iter << " - Error: " << error << std::endl;
        }
    }

    std::cout << "Jacobi terminated. Iterations: " << iter << ", Final Error: " << error << std::endl;
}

double SerialSolver::compute_analytical_error(std::function<double(double, double)> exact_sol) {
    double total_error_sum = 0.0;

    // Cicla su TUTTA la griglia globale
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double x = i * h;
            double y = j * h;

            double diff = U(i, j) - exact_sol(x, y);
            total_error_sum += diff * diff;
        }
    }

    // Applica la formula del PDF direttamente alla fine
    return std::sqrt(h * total_error_sum);
}
#pragma once
#include "Matrix.hpp"
#include <functional>

/**
 * @brief SerialSolver class encapsulates the logic for solving the 2D Poisson equation in serial.
 */
class SerialSolver {
private:
    int n;
    double h;       
    
    Matrix U; // Grid at step k
    Matrix U_new; // Grid at step k+1
    
    std::function<double(double, double)> f;
    std::function<double(double, double)> bc;

public:
    /**
     * @brief Constructor for the SerialSolver class.
     * @param grid_size The size of the grid.
     * @param forcing_term The forcing term function.
     * @param boundary_term The boundary condition function.
     */
    SerialSolver(int grid_size, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term);

    /**
     * @brief Solves the 2D Poisson equation in serial.
     * @param max_iterations The maximum number of iterations.
     * @param tolerance The tolerance for convergence.
     */
    void solve(int max_iterations, double tolerance);

    /**
     * @brief Computes the analytical error of the numerical solution.
     * @param exact_sol The exact solution function.
     * @return The computed error.
     */
    double compute_analytical_error(std::function<double(double, double)> exact_sol) const ;
};
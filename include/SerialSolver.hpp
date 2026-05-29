#pragma once
#include "Matrix.hpp"
#include <functional>

class SerialSolver {
private:
    int n;
    double h;       
    
    Matrix U; // Grid at step k
    Matrix U_new; // Grid at step k+1
    
    std::function<double(double, double)> f;
    std::function<double(double, double)> bc;

public:
    // Constructor
    SerialSolver(int grid_size, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term);

    void solve(int max_iterations, double tolerance);

    double compute_analytical_error(std::function<double(double, double)> exact_sol);
};
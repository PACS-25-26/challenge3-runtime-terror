#pragma once
#include "Matrix.hpp"
#include "MpiDomain.hpp"
#include <functional>

class ParallelSolver {
private:
    MpiDomain domain;
    double h;
    
    Matrix U;       
    Matrix U_new;   
    
    std::function<double(double, double)> f;

public:
    // Constructor takes the MPI domain configuration
    ParallelSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term);
    
    // Parallel solve method
    void solve(int max_iter, double tol);
};
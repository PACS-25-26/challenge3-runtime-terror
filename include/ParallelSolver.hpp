#pragma once
#include "Matrix.hpp"
#include "MpiDomain.hpp"
#include <functional>
#include<string>

class ParallelSolver {
private:
    MpiDomain domain;
    double h;
    
    Matrix U;       
    Matrix U_new;   
    
    std::function<double(double, double)> f;
    std::function<double(double, double)> bc;

public:
    // Constructor takes the MPI domain configuration
    ParallelSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term);
    
    // Parallel solve method
    void solve(int max_iter, double tol);

    double compute_analytical_error(std::function<double(double, double)> exact_sol);

    void export_vtk(const std::string& filename); 
};
#pragma once
#include "Matrix.hpp"
#include "MpiDomain.hpp"
#include "VTKExport.hpp"
#include <functional>
#include<string>

/**
 * @brief ParallelSolver class encapsulates the logic for solving the 2D Poisson equation in parallel using MPI.
 */
class ParallelSolver {
private:
    MpiDomain domain;
    double h;
    
    Matrix U;       
    Matrix U_new;   
    
    std::function<double(double, double)> f;
    std::function<double(double, double)> bc;

public:
    /**
     * @brief Constructor for the ParallelSolver class.
     * * Initializes the solver with the given MPI domain, forcing term, and boundary condition functions
     */
    ParallelSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term);
    
    /**
     * @brief Solves the 2D Poisson equation in parallel using MPI.
     * @param max_iter Maximum number of iterations.
     * @param tol Tolerance for convergence.
     */
    void solve(int max_iter, double tol);

    /**
     * @brief Computes the analytical error of the numerical solution.
     * @param exact_sol The exact solution function.
     * @return The computed error.
     */
    double compute_analytical_error(std::function<double(double, double)> exact_sol) const;

    /**
     * @brief Exports the solution to a VTK file.
     * @param filename The name of the output file.
     */
    void export_vtk(const std::string& filename){
        VTKExport exporter(domain, h);
        exporter.export_vtk(U, filename);
    }
};
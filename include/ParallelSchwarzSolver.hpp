#pragma once
#include "Matrix.hpp"
#include "MpiDomain.hpp"
#include "VTKExport.hpp"
#include <functional>
#include<string>

/**
 * @brief parallel Block Jacobi solver for the 2D Laplace equation
 * 
 */
class ParallelSchwarzSolver {
private:
    MpiDomain domain;
    double h;

    Matrix U;       
    Matrix U_new;   
    
    std::function<double(double, double)> f;
    std::function<double(double, double)> bc;

public:
    /**
     * @brief constructor for the ParallelSchwarzSolver class
     * * Initializes the solver with the given MPI domain, forcing term, and boundary condition functions
     */
    ParallelSchwarzSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term);
    

    /**
    * @brief Implements parallel block Jacobi iterations: the domain is decomposed into subdomains (blocks) and each block is solved iteratively while exchanging BC info with the neighbors 
    * 
    * @param max_global_iter maximum number of block Jacobi iterations (global iterations)
    * @param max_local_iter maximum number of local iterations
    * @param tol tolerance for convergence
    */
    void solve(int max_global_iter, int max_local_iter, double tol);


    /**
    * @brief Computes the L2 error between the numerical solution and the exact solution
    * 
    * @param exact_sol  exact solution function
    * @return L2 error
    */
    double compute_analytical_error(std::function<double(double, double)> exact_sol) const;


    /**
     * @brief Exports the solution to a VTK file
     * @param filename name of the output file
     */
    void export_vtk(const std::string& filename){
        VTKExport exporter(domain, h);
        exporter.export_vtk(U, filename);
    } 
};
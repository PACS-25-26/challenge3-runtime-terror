#include "ParallelSchwarzSolver.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h> 
#include <omp.h> 


ParallelSchwarzSolver::ParallelSchwarzSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term, std::function<double(double, double)> boundary_term)
    : domain(mpi_dom), 
      h(1.0 / (domain.global_n - 1.0)), 
      // Allocate using allocated_rows (local rows + 2 ghost cells)
      U(domain.allocated_rows, domain.global_n),       
      U_new(domain.allocated_rows, domain.global_n),   
      f(forcing_term),
      bc(boundary_term) 
{
    for(int i = 1; i <= domain.local_rows; ++i){
        int global_i = domain.start_row + (i - 1);
        double x = global_i * h;

        for(int j = 0; j < domain.global_n; ++j){
            double y = j * h;

            if(global_i == 0 || global_i == domain.global_n - 1 || j == 0 || j == domain.global_n - 1){
                U(i, j) = bc(x, y);
                U_new(i, j) = bc(x, y);
            }
        }
    }
}


void ParallelSchwarzSolver::solve(int max_global_iter, int max_local_iter, double tol){
    double global_error = tol + 1.0;
    int global_iter = 0; 

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int top_neighbor = (domain.start_row > 0) ? rank - 1 : MPI_PROC_NULL;
    int bottom_neighbor = (domain.start_row + domain.local_rows < domain.global_n) ? rank + 1 : MPI_PROC_NULL;

    if (rank == 0){
        std::cout << "Starting MPI Block Jacobi iterations..." << std::endl;
    }

    while(global_error > tol && global_iter < max_global_iter){
        MPI_Sendrecv(&U(1, 0), domain.global_n, MPI_DOUBLE, top_neighbor, 0,
                    &U(domain.local_rows +1, 0), domain.global_n, MPI_DOUBLE, bottom_neighbor, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&U(domain.local_rows, 0), domain.global_n, MPI_DOUBLE, bottom_neighbor, 1,
                    &U(0, 0), domain.global_n, MPI_DOUBLE, top_neighbor, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if(top_neighbor != MPI_PROC_NULL){
            for(int j = 0; j < domain.global_n; ++j){
                U_new(0, j) = U(0, j); }
        }
        if(bottom_neighbor != MPI_PROC_NULL){
            for(int j = 0; j < domain.global_n; ++j){
                U_new(domain.local_rows + 1, j) = U(domain.local_rows + 1, j); }
        }

        // local solver
        double block_diff_sum = 0.0;
        double local_residual = tol + 1.0;
        int local_iter = 0;

        Matrix U_old_block = U; // local copy for each block

        while(local_residual > tol && local_iter < max_local_iter){
            local_residual = 0.0;

            int i_start = (domain.start_row == 0) ? 2 : 1; 
            int i_end = (domain.start_row + domain.local_rows == domain.global_n) ? domain.local_rows - 1 : domain.local_rows;

            #pragma omp parallel for reduction(+:local_residual)
            for(int i = i_start; i <= i_end; ++i){
                for(int j = 1; j < domain.global_n - 1; ++j){
                    int global_i = domain.start_row + (i-1);
                    double x = global_i * h;
                    double y = j * h;

                    U_new(i, j) = 0.25 * (U(i - 1, j) + U(i + 1, j) + U(i, j - 1) + U(i, j + 1) + (h * h * f(x, y)));

                    double diff = U_new(i, j) - U(i, j);
                    local_residual += diff * diff;
                }
            }

            local_residual = std::sqrt(h * local_residual);
            U.swap(U_new);
            local_iter++;
    }

        // global error
        for(int i = 1; i <= domain.local_rows; ++i){
            for(int j = 1; j < domain.global_n -1; ++j){
                double diff = U(i, j) - U_old_block(i, j);
                block_diff_sum += diff * diff;
            }
        }
        
        double global_diff_sum = 0.0;
        MPI_Allreduce(&block_diff_sum, &global_diff_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        global_error = std::sqrt(h * global_diff_sum);
        global_iter++;

        if(rank == 0 && global_iter % 100 == 0){
            std::cout << "Global Block Jacobi Iteration " << global_iter << " - Global Error: " << global_error << std::endl;
        }
    }

    if (rank == 0){
        std::cout << "MPI Block Jacobi terminated. Global Iterations: " << global_iter 
                  << ", Final Global Error: " << global_error << std::endl;
    }
}


double ParallelSchwarzSolver::compute_analytical_error(std::function<double(double, double)> exact_sol) const{
    double local_error_sum = 0.0;

    #pragma omp parallel for reduction(+:local_error_sum)
    for (int i = 1; i <= domain.local_rows; ++i) {
        int global_i = domain.start_row + (i - 1);

        for (int j = 0; j < domain.global_n; ++j) {
            double x = global_i * h;
            double y = j * h;

            double diff = U(i, j) - exact_sol(x, y);
            local_error_sum += diff * diff;
        }
    }

    double global_error_sum = 0.0;
    
    // Reduce all local sums to the master process 
    // We use MPI_Reduce instead of Allreduce because only Rank 0 needs to print it
    MPI_Reduce(&local_error_sum, &global_error_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // The L2 norm over a 2D grid is scaled by the area element (h * h)
    // Error = sqrt( sum(diff^2) * h^2 ) = h * sqrt(sum(diff^2))
    return std::sqrt(h *global_error_sum);
}


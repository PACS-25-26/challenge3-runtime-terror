#include "ParallelSolver.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h> // We need MPI here for the communication
#include <omp.h> // We need OpenMP for parallel loops if we choose to use it in the future

// Constructor
ParallelSolver::ParallelSolver(MpiDomain mpi_dom, std::function<double(double, double)> forcing_term)
    : domain(mpi_dom), 
      h(1.0 / (domain.global_n - 1.0)), 
      // Allocate using allocated_rows (local rows + 2 ghost cells)
      U(domain.allocated_rows, domain.global_n),       
      U_new(domain.allocated_rows, domain.global_n),   
      f(forcing_term) 
{
}

void ParallelSolver::solve(int max_iter, double tol) {
    double global_error = tol + 1.0; 
    int iter = 0;

    // Identify neighbors. If a rank is at the true boundary, it has no neighbor (MPI_PROC_NULL)
    int top_neighbor = (domain.rank == 0) ? MPI_PROC_NULL : domain.rank - 1;
    int bottom_neighbor = (domain.rank == domain.size - 1) ? MPI_PROC_NULL : domain.rank + 1;

    if (domain.rank == 0) {
        std::cout << "Starting MPI Jacobi iterations..." << std::endl;
    }

    while (global_error > tol && iter < max_iter) {
        double local_error_sum = 0.0;
        
        // Send top real row (1) to top_neighbor
        // Receive into bottom ghost cell (allocated_rows - 1) from bottom_neighbor
        MPI_Sendrecv(
            &U(1, 0), domain.global_n, MPI_DOUBLE, top_neighbor, 0,
            &U(domain.allocated_rows - 1, 0), domain.global_n, MPI_DOUBLE, bottom_neighbor, 0,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );

        // Send bottom real row (allocated_rows - 2) to bottom_neighbor
        // Receive into top ghost cell (0) from top_neighbor
        MPI_Sendrecv(
            &U(domain.allocated_rows - 2, 0), domain.global_n, MPI_DOUBLE, bottom_neighbor, 1,
            &U(0, 0), domain.global_n, MPI_DOUBLE, top_neighbor, 1,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );
        
        #pragma omp parallel for reduction(+:local_error_sum) 
        // Loop only over the actual rows assigned to this rank (from 1 to local_rows)
        for (int i = 1; i <= domain.local_rows; ++i) {
            
            // Compute the global row index to check for physical boundaries and coordinates
            int global_i = domain.start_row + (i - 1);
            
            // If this row is the physical top (0) or bottom (n-1) of the whole grid, skip it (Dirichlet BC already set)
            if (global_i == 0 || global_i == domain.global_n - 1) {
                continue; 
            }

            for (int j = 1; j < domain.global_n - 1; ++j) {
                double x = global_i * h;
                double y = j * h; 

                // Update formula
                U_new(i, j) = 0.25 * (U(i - 1, j) + U(i + 1, j) + U(i, j - 1) + U(i, j + 1) + (h * h * f(x, y)));

                double diff = U_new(i, j) - U(i, j);
                local_error_sum += diff * diff;
            }
        }
        
        double global_error_sum = 0.0;
        // Sum all local_error_sum into global_error_sum across all ranks
        MPI_Allreduce(&local_error_sum, &global_error_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // Compute the final L2 norm error
        global_error = std::sqrt(h * global_error_sum);

        U.swap(U_new); // O(1) pointer swap
        
        iter++;
        
        if (domain.rank == 0 && iter % 1000 == 0) {
            std::cout << "Iteration " << iter << " - Error: " << global_error << std::endl;
        }
    }

    if (domain.rank == 0) {
        std::cout << "MPI Jacobi terminated. Iterations: " << iter 
                  << ", Final Error: " << global_error << std::endl;
    }
}

double ParallelSolver::compute_analytical_error(std::function<double(double, double)> exact_sol) {
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


void ParallelSolver::export_vtk(const std::string& filename) {
    int n = domain.global_n;

    std::vector<int> counts(domain.size);
    std::vector<int> displs(domain.size);

    int base_rows = n / domain.size;
    int remainder = n % domain.size;

    for(int i = 0; i < domain.size; ++i){
        int rows_rank_i = (i < remainder) ? base_rows + 1 : base_rows;
        counts[i] = rows_rank_i * n;

        int start_row_i = (i < remainder) ? i * rows_rank_i : remainder * (base_rows + 1) + (i - remainder) * base_rows;
        displs[i] = start_row_i * n;
    }

    std::vector<double> global_U;
    if (domain.rank == 0) {
        global_U.resize(n * n, 0.0);
    }
    
    MPI_Gatherv(&U(1, 0), domain.local_rows * n, MPI_DOUBLE, 
                global_U.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                0, MPI_COMM_WORLD);

    if (domain.rank == 0) {
        std::ofstream vtk_file(filename);
        if(!vtk_file){
            std::cerr << "Error writing VTK file" << std::endl;
            return;
        }
        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << "Laplace Solution\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_POINTS\n";
        vtk_file << "DIMENSIONS " << n << " " << n << " 1\n";
        vtk_file << "ORIGIN 0 0 0\n";
        vtk_file << "SPACING " << h << " " << h << " 0\n";
        vtk_file << "POINT_DATA " << n * n << "\n";
        vtk_file << "SCALARS U double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";

        for(int i = 0; i < n * n; ++i){
            vtk_file << global_U[i] << "\n";
        }
        vtk_file.close();
        std::cout << "VTK file '" << filename << "' written successfully." << std::endl;
    }
}
#include "VTKExport.hpp"
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>

VTKExport::VTKExport(MpiDomain mpi_dom, double spacing) : domain(mpi_dom), h(spacing) {}

/**
 * @brief Exports the solution to VTK format. Each process sends local data to the root process, which writes the complete solution in the VTK file.
 * 
 * @param U solution matrix
 * @param filename name of the output VTK file
 */
void VTKExport::export_vtk(const Matrix& U, const std::string& filename) {
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
#pragma once
#include "Matrix.hpp"
#include "MpiDomain.hpp"
#include <string>

/**
 * @brief Exports the solution of the PDE to VTK file for visualization in Paraview
 * 
 */
class VTKExport{
    private:
        MpiDomain domain;
        double h;
    public:
        VTKExport(MpiDomain mpi_dom, double spacing);

        void export_vtk(const Matrix& U, const std::string& filename);
};
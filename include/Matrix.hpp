#pragma once
#include <vector>

/**
 * @brief A simple 2D dense matrix class flattened into a 1D vector.
 * * This class provides basic matrix operations, element access via (i, j) 
 * coordinates (using row-major order), and a raw pointer access which is 
 * particularly useful for MPI communications.
 */
class Matrix {
private:
    int rows;
    int cols;
    std::vector<double> data;

public:
    /**
     * @brief Constructor for the Matrix class.
     * * Allocates memory for the matrix and initializes all elements to zero.
     * @param r Number of rows.
     * @param c Number of columns.
     */
    Matrix(int r, int c) : rows(r), cols(c), data(r * c, 0.0) {}

    /**
     * @brief Get the number of rows in the matrix.
     * @return Number of rows.
     */
    int getRows() const { return rows; }

    /**
     * @brief Get the number of columns in the matrix.
     * @return Number of columns.
     */
    int getCols() const { return cols; }

    /**
     * @brief Access element at position (i, j) for writing.
     * @param i Row index.
     * @param j Column index.
     * @return Reference to the element at (i, j).
     */
    double& operator()(int i, int j) {
        return data[i * cols + j];
    }

    /**
     * @brief Access element at position (i, j) for reading.
     * @param i Row index.
     * @param j Column index.
     * @return Const reference to the element at (i, j).
     */
    const double& operator()(int i, int j) const {
        return data[i * cols + j];
    }

    /**
     * @brief Returns a pointer to the underlying data array.
     * * This is particularly useful for MPI communications.
     * @return Pointer to the first element of the data array.
     */
    double* data_ptr() {
        return data.data();
    }

        /**
        * @brief Returns a const pointer to the underlying data array.
        * @return Const pointer to the first element of the data array.
        */
    const double* data_ptr() const {
        return data.data();
    }

    /**
     * @brief Swaps the contents of this matrix with another matrix.
     * @param other Matrix to swap with.
     */
    void swap(Matrix& other) {
        std::swap(rows, other.rows);
        std::swap(cols, other.cols);
        std::swap(data, other.data);
    }
};
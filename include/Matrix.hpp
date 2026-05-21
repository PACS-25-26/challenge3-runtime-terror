#pragma once
#include <vector>

class Matrix {
private:
    int rows;
    int cols;
    std::vector<double> data;

public:
    // Constructor: allocate memory for the matrix and initialize to zero
    Matrix(int r, int c) : rows(r), cols(c), data(r * c, 0.0) {}

    // Getters
    int getRows() const { return rows; }
    int getCols() const { return cols; }

    // Write access
    double& operator()(int i, int j) {
        return data[i * cols + j];
    }

    // Read only access
    const double& operator()(int i, int j) const {
        return data[i * cols + j];
    }

    // Returns pointer to base vector
    double* data_ptr() {
        return data.data();
    }

    void swap(Matrix& other) {
        std::swap(rows, other.rows);
        std::swap(cols, other.cols);
        std::swap(data, other.data);
    }
};
#pragma once
#include <functional>

/**
 * @brief Struct to hold the problem data for PDE solvers, including forcing term, exact solution and BCs
 * 
 */
struct ProblemData{
    /**
     * @brief Forcing term f(x, y) in the PDE
     */
    std::function<double(double, double)> f;
    /**
     * @brief Exact solution u(x, y) of the PDE
     */
    std::function<double(double, double)> exact_solution;
    /**
     * @brief Boundary condition g(x, y) on the domain boundary
     */
    std::function<double(double, double)> g;
};

/**
 * @brief Namespace containing standard test cases for the Laplace equation.
 * * Provides predefined test cases with homogeneous and inhomogeneous Dirichlet boundary conditions.
 */
namespace TestCases {
    /**
     * @brief Test case from the reference paper with homogeneous Dirichlet BCs.
     * * The exact solution is u(x,y) = sin(2*pi*x)*sin(2*pi*y).
     * @return ProblemData object populated with the respective functions.
     */
    ProblemData Sine(); 
    /**
     * @brief Test case with inhomogeneous Dirichlet BCs.
     * * The exact solution is u(x,y) = exp(x) * cos(y).
     * @return ProblemData object populated with the respective functions.
     */
    ProblemData ExpCos(); 
} // end of namespace TestCases
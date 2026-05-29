#pragma once
#include <functional>

/**
 * @brief struct to hold the problem data for PDE solvers, including forcing term, exact solution and BCs
 * 
 */
struct ProblemData{
    // Forcing term f(x, y)
    std::function<double(double, double)> f;
    // Exact solution u(x, y)
    std::function<double(double, double)> exact_solution;
    // Boundary condition g(x, y)
    std::function<double(double, double)> g;
};

/**
 * @brief Test cases for homogenuous and inhomogenuous Dirichlet BCs
 * 
 */
namespace TestCases {
    ProblemData Sine(); // paper test case with u(x,y) = sin(2*pi*x)*sin(2*pi*y)
    ProblemData ExpCos(); // test case with u(x,y) = exp(x) * cos(y);
} // end of namespace TestCases
#pragma once
#include <functional>

struct ProblemData{
    // Forcing term f(x, y)
    std::function<double(double, double)> f;
    // Exact solution u(x, y) (useful later for computing L2 error)
    std::function<double(double, double)> exact_solution;
    // Boundary condition g(x, y)
    std::function<double(double, double)> g;
};

namespace TestCases {
    ProblemData Sine(); // paper test case with u(x,y) = sin(2*pi*x)*sin(2*pi*y)
    ProblemData ExpCos(); // test case with u(x,y) = exp(x) * cos(y);
} // end of namespace TestCases
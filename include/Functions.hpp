#pragma once

namespace EquationData {
    // Forcing term f(x, y)
    double forcing_term(double x, double y);

    // Exact solution u(x, y) (useful later for computing L2 error)
    double exact_solution(double x, double y);
}
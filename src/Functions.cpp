#include "Functions.hpp"
#include <cmath>

namespace EquationData {

    double forcing_term(double x, double y) {
        return 8.0 * M_PI * M_PI * std::sin(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y);
    }

    double exact_solution(double x, double y) {
        return std::sin(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y);
    }

}
#include "Functions.hpp"
#include <cmath>

namespace TestCases{
    ProblemData Sine(){
        return{
            // f(x,y)
            [](double x, double y){ return 2.0 * M_PI * 2.0 * M_PI * std::sin(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y); },
            // exact(x,y)
            [](double x, double y){ return std::sin(2.0 * M_PI * x) * std::sin(2.0 * M_PI * y); },
            // g(x,y)
            [](double x, double y){ return 0.0 * x + 0.0 * y; }
        };
    }

    ProblemData ExpCos(){
        return{
            // f(x,y)
            [](double x, double y){return 0.0 * x + 0.0 * y;},
            // exact(x,y)
            [](double x, double y){ return std::exp(x) * std::cos(y); },
            // g(x,y)
            [](double x, double y){ return std::exp(x) * std::cos(y); }     
        };
    }
} // end of namespace TestCases

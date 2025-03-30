#include "abc.h"

extern "C" {
    double* solve_c(double* numbers, double* chebyshev, int N_, int M_) {
        return solve(numbers, chebyshev, N_, M_);
    }
}
#ifndef TEST_H
#define TEST_H
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <limits>
#include <random>
#include <chrono>


double* solve(double* numbers, double* chebyshev, int N_, int M_);

extern "C" {
    double* solve_c(double* numbers, double* chebyshev, int N_, int M_);
}
#endif
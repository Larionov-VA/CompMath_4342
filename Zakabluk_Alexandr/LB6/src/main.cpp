#include <iostream>
#include <cmath>
#include "methods.hpp"

#ifndef X0_START
    #define X0_START 1.0 // Начальное приближение x0
#endif

int main() {
    int iterationsCount;
    std::cout << "Eps\t\t | Root\t\t\t | Iterations\n";
    std::cout << "--------------------------------------------------\n";
    for (double epsilon = 0.1; epsilon >= 0.000001; epsilon /= 10) {
        std::cout << epsilon << "\t\t | "
                  << ITER(X0_START, epsilon, iterationsCount) << "\t\t | "
                  << iterationsCount << '\n';
    }
    return 0;
}
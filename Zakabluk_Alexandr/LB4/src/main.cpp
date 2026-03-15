#include <iostream>
#include <cmath>
#include "methods.hpp"

#ifndef LEFT
    #define LEFT 0.5 // левая граница
#endif
#ifndef RIGHT
    #define RIGHT 1.5 // правая граница
#endif

int main() {
    int iterationsCount;
    std::cout << "Eps\t\t | Root\t\t\t | Iterations\n";
    std::cout << "--------------------------------------------------\n";
    for (double epsilon = 0.1; epsilon >= 0.000001; epsilon /= 10) {
        std::cout << epsilon << "\t\t | "
                  << HORDA(LEFT, RIGHT, epsilon, iterationsCount) << "\t\t | "
                  << iterationsCount << '\n';
    }
    return 0;
}
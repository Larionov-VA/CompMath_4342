#include <iostream>
#include <cmath>
#include <iomanip>
#include "methods.hpp"

#ifndef LEFT
#define LEFT 2.8
#endif
#ifndef RIGHT
#define RIGHT 3.2
#endif

int main() {
    int iterationsCount;
    std::cout << "Epsilon\t\t | Root (x)\t | Iterations\n";
    std::cout << "--------------------------------------------------\n";
    
    for (double epsilon = 0.1; epsilon >= 0.000001; epsilon /= 10) {
        double root = BISECT(LEFT, RIGHT, epsilon, iterationsCount);
        std::cout << std::fixed << std::setprecision(6) 
                  << epsilon << "\t | " 
                  << root << "\t | " 
                  << iterationsCount << '\n';
    }
    return EXIT_SUCCESS;
}
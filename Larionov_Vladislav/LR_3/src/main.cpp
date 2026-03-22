#include <iostream>
#include <cmath>
#include "../../methods.hpp"

#ifndef LEFT
    // минимум для монотонности -- 0.84472
    #define LEFT 1.5L
#endif
#ifndef RIGHT
    // максимум -- 9223372036854775807 (максимальный double)
    #define RIGHT 2.L
#endif
#ifndef EPSILON
    int main() {
        int iterationsCount;
        for (double epsilon = 0.1; epsilon >= 0.000001; epsilon /= 10) {
            std::cout <<
                epsilon << "\t | " <<
                NEWTON(RIGHT, epsilon, iterationsCount) << "\t | " <<
                iterationsCount << '\n';
        }
        return EXIT_SUCCESS;
    }
#else
    int main() {
        int iterationsCount;
        std::cout << NEWTON(RIGHT, EPSILON, iterationsCount) << "\t | " << iterationsCount << '\n';
        return EXIT_SUCCESS;
    }
#endif
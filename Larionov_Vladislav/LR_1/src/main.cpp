#include <iostream>
#include <cmath>
#include "../../methods.hpp"

#ifndef LEFT
    #define LEFT 1.5L // минимум для монотонности -- 0.84472
#endif
#ifndef RIGHT
    #define RIGHT 2.L // максимум -- 9223372036854775807 (максимальный double)
#endif

int main() {
    int iterationsCount;
    for (double epsilon = 0.1; epsilon >= 0.000001; epsilon /= 10) {
        std::cout <<
            epsilon << "\t | " <<
            BISECT(LEFT, RIGHT, epsilon, iterationsCount) << "\t | " <<
            iterationsCount << '\n';
    }
    return EXIT_SUCCESS;
}
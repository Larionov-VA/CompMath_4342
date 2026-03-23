/*
 * ЛР -- 4. Вариант -- 8:
 * f(x) = 2x^2 - x^4 - 1 - ln(x)
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "methods.h"

double F(double x) {
    return 2 * pow(x, 2) - pow(x, 4) - 1 - log(x);
}

int main() {
    int n;
    double x;
    double eps;
    double delta;

    double a = 0.85;
    double b = 1.15;

    for (int di=1; di<7; di++) {
        delta = pow(10, -di);
        printf("For Delta = %lf:\n", delta);

        for (int i=1; i<7; i++) {
            eps = pow(10, -i);
            x = HORDA(a, b, eps, n);

            printf("\tFor epsilon = %lf: x = %lf, iter = %d\n", eps, Round(x, delta), n);
        }
    }

    return 0;
}
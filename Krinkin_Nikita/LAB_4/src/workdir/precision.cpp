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

int main(int argc, char *argv[]) {
    int n;
    double x;
    double eps;
    double delta = pow(10, -4);

    FILE *out_file = fopen(argv[1], "w");

    double a = 0.85;
    double b = 1.15;

    for (int i=1; i<7; i++) {
        eps = pow(10, -i);
        x = HORDA(a, b, eps, n);
        fprintf(out_file, "%d %d\n", i, n);
    }
    fclose(out_file);

    return 0;
}
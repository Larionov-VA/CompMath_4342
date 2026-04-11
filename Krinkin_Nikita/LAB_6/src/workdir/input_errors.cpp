/*
 * ЛР -- 6. Вариант -- 8:
 * f(x) = 2x^2 - x^4 - 1 - ln(x)
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "methods.h"

using namespace std;

double origF(double x) {
    return 2 * pow(x, 2) - pow(x, 4) - 1 - log(x);
}

double F(double x) {
    return pow(exp(1.0), 2 * pow(x, 2) - pow(x, 4) - 1);
}

int main(int argc, char *argv[]) {
    int n;
    double x;
    double eps = pow(4, -10);

    FILE *out_file = fopen(argv[1], "w");

    for (double x0=0; x0 < 2; x0 += 0.001) {
        x = ITER(x0, eps, n);
        fprintf(out_file, "%lf %d\n", x0, n);
    }

    return 0;
}
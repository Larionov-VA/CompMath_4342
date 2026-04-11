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
    double eps;
    double delta = pow(10, -4);

    FILE *out_file = fopen(argv[1], "w");

    double x0 = 1.05;

    for (int i=1; i<=12; i++) {
        eps = pow(10, -i);
        x = ITER(x0, eps, n);
        fprintf(out_file, "%d %d\n", i, n);
    }
    fclose(out_file);

    return 0;
}

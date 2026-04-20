/*
 * ЛР -- 5. Вариант -- 8:
 * f(x) = 2x^2 - x^4 - 1 - ln(x)
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#define __NEWTON

#include "methods.h"

using namespace std;

double F(double x) {
    return 2 * pow(x, 2) - pow(x, 4) - 1 - log(x);
}

double F1(double x) {
    return 4 * x - 4 * pow(x, 3) - pow(x, -1);
}

int main(int argc, char *argv[]) {
    int n;
    double x;
    double eps = pow(8, -10);

    FILE *out_file = fopen(argv[1], "w");

    for (double x0=0.2; x0 < 2; x0 += 0.001) {
        x = NEWTON(x0, eps, n);
        fprintf(out_file, "%lf %d\n", x0, n);
    }

    return 0;
}
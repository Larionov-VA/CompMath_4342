/*
 * ЛР -- 3. Вариант -- 8:
 * f(x) = 2x^2 - x^4 - 1 - ln(x)
 */

#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

#include "methods.h"

// using namespace std;


double F(double x) {
    return 2 * pow(x, 2) - pow(x, 4) - 1 - log(x);
}

int main(int argc, char* argv[]) {
    int n;
    double x;
    double eps;
    double delta = pow(10, -4);
    
    FILE *out_file = fopen(argv[1], "w");

    double a = 0.13453;
    double b = 1.4432341923;

    for (int i=1; i<7; i++) {
        eps = pow(10, -i);
        x = BISECT(a, b, eps, n);

        printf("\tFor epsilon = %lf: x = %lf, iter = %d\n", eps, Round(x, delta), n);

        fprintf(out_file, "%d %d\n", i, n);
    }

    fclose(out_file);

    return 0;
}
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

int main() {
    double x;
    double eps;
    double delta;
    
    int it;
    cout << "Enter epsilon presicion (decimal points): ";
    cin >> it; if (it <= 0) it = 4;
    eps = pow(10, -it);

    cout << "Enter delta presicion (decimal points): ";
    cin >> it; if (it <= 0) it = 4;
    delta = pow(10, -it);

    double x0 = 1.05;
    int n = 0;
    x = NEWTON(x0, eps, n);

    printf("\tFor epsilon = %lf and x0 = %lf: x = %lf with delta = %lf, iter = %d\n", eps, x0, Round(x, delta), delta, n);

    return 0;
}
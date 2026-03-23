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

using namespace std;


double F(double x) {
    return 2 * pow(x, 2) - pow(x, 4) - 1 - log(x);
}

int main() {
    double x;
    double eps;
    double delta;
    
    int it;
    cout << "Enter epsilon presicion (decimal points): ";
    cin >> it; if (it <= 0) it = 1;
    eps = pow(10, -it);

    cout << "Enter delta presicion (decimal points): ";
    cin >> it; if (it <= 0) it = 1;
    delta = pow(10, -it);

    double a = 0.1;
    double b = 1.5;

    int n = 0;
    x = HORDA(a, b, eps, n);

    printf("\tFor epsilon = %lf: x = %lf with delta = %lf, iter = %d\n", eps, Round(x, delta), delta, n);

    return 0;
}
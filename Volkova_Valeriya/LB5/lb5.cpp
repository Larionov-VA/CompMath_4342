#include <iostream>
#include <cmath>
#include <iomanip>

#define __NEWTON
#include "methods.h"

using namespace std;

double global_delta = 0.0;

double F(double x)
{
    double y = acos(x * x) - x;
    if (global_delta > 0.0) return Round(y, global_delta);
    return y;
}

double F1(double x)
{
    double y = -(2.0 * x) / sqrt(1.0 - x * x * x * x) - 1.0;
    if (global_delta > 0.0) return Round(y, global_delta);
    return y;
}

int main()
{
    double x0 = 0.9;
    int N;

    cout << fixed << setprecision(10);

    cout << "Зависимость N от Eps" << endl;
    for (double eps = 0.1; eps >= 1e-6; eps /= 10.0)
    {
        double root = NEWTON(x0, eps, N);
        cout << "Eps: " << eps
             << "\tN: " << N
             << "\tRoot: " << root
             << "\tF(root): " << F(root) << endl;
    }

    cout << endl << "Влияние Delta" << endl;
    global_delta = 1e-7;
    cout << "Delta: " << global_delta << endl;
    for (double eps = 0.1; eps >= 1e-6; eps /= 10.0)
    {
        double root = NEWTON(x0, eps, N);
        cout << "Eps: " << eps
                << "\tN: " << N
                << "\tRoot: " << root
                << "\tF(root): " << F(root) << endl;
    }


    return 0;
}

#include <iostream>
#include <cmath>
#include <iomanip>

#include "methods.h"

using namespace std;

double global_delta = 0.0;

double Equation(double x)
{
    return acos(x * x) - x;
}

double F(double x)
{
    double y = sqrt(cos(x));
    if (global_delta > 0.0) return Round(y, global_delta);
    return y;
}

double F1(double x)
{
    double y = -sin(x) / (2.0 * sqrt(cos(x)));
    if (global_delta > 0.0) return Round(y, global_delta);
    return y;
}

int main()
{
    double x0 = 0.8;
    int N;

    cout << fixed << setprecision(10);

    cout << "Зависимость N от Eps" << endl;
    global_delta = 0.0;
    for (double eps = 0.1; eps >= 1e-6; eps /= 10.0)
    {
        double root = ITER(x0, eps, N);
        cout << "Eps: " << eps
             << "\tN: " << N
             << "\tRoot: " << root
             << "\tf(root): " << Equation(root) << endl;
    }

    cout << endl << "Влияние Delta" << endl;
    double deltas[] = {1e-1, 1e-3, 1e-5, 1e-7};

    for (double delta : deltas)
    {
        global_delta = delta;
        cout << "Delta: " << global_delta << endl;

        for (double eps = 0.1; eps >= 1e-6; eps /= 10.0)
        {
            double root = ITER(x0, eps, N);
            cout << "Eps: " << eps
                 << "\tN: " << N
                 << "\tRoot: " << root
                 << "\tf(root): " << Equation(root) << endl;
        }

        cout << endl;
    }

    return 0;
}

#include <iostream>
#include <cmath>
#include "methods.h"

using namespace std;


double global_delta = 0.0;

double F(double x) {
    double y = acos(x*x) - x;
    if (global_delta > 0) return Round(y, global_delta);
    return y;
}


int main() {
    double left = 0.0, right = 1.0;
    int N;

    cout << "Зависимость N от Eps" << endl;
    for (double eps = 0.1; eps >= 1e-6; eps /= 10) {
        double root = BISECT(left, right, eps, N);
        cout << "Eps:" << eps << "\tN:" << N << "\tRoot:" << root << endl;
    }

    cout << "Влияние Delta" << endl;
    for (double delta = 0.1; delta >= 1e-7; delta /= 10) {
        global_delta = delta;
        cout << "Delta: " << delta << endl;
        for (double eps = 0.1; eps >= 1e-6; eps /= 10) {
            double root = BISECT(left, right, eps, N);
            cout << "Eps:" << eps << "\tN:" << N << "\tRoot:" << root << endl;
        }
    }

    
    return 0;
}
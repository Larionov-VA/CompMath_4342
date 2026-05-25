#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double Delta = 0.0;

double Round(double X, double Delta) {
    if (Delta <= 1e-12) return X;
    return Delta * round(X / Delta);
}

double F(double x) {
    double y = cos(x * x) / 10.0;
    if (Delta > 1e-12) return Round(y, Delta);
    return y;
}

double F1(double x) {
    return -0.2 * x * sin(x * x);
}

double ITER(double x0, double eps, int &iter, int max_iter = 1000) {
    double x = x0;
    double x_prev;
    iter = 0;
    do {
        x_prev = x;
        x = F(x_prev);
        iter++;
        if (iter >= max_iter) {
            cerr << "Превышен лимит итераций (возможно, зацикливание) при Delta = " << Delta << endl;
            break;
        }
    } while (fabs(x - x_prev) >= eps);
    return x;
}

int main() {
    double x0 = 1.0;
    int iter;

    cout << "=== Зависимость числа итераций от eps ===" << endl;
    cout << fixed << setprecision(12);
    for (double eps = 0.1; eps >= 1e-6; eps /= 10) {
        double root = ITER(x0, eps, iter);
        cout << eps << "\t" << root << "\t" << iter << endl;
    }

    cout << "\n=== Влияние округления (eps=1e-6) ===" << endl;
    double eps_fixed = 1e-6;
    double Deltas[] = {0.0, 0.001, 0.01, 0.1, 1.0};
    for (double d : Deltas) {
        Delta = d;
        double root = ITER(x0, eps_fixed, iter);
        cout << d << "\t" << root << "\t" << iter << endl;
    }
    return 0;
}
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

double Delta = 0.0;

double f(double x) {
    return cos(x * x) - 10.0 * x;
}

double df(double x) {
    return -2.0 * x * sin(x * x) - 10.0;
}

double Round(double X, double Delta) {
    if (Delta <= 1e-12) return X;
    return Delta * round(X / Delta);
}

double newton(double x0, double eps, int &iter, bool use_round = false, int max_iter = 100) {
    double x = x0;
    double x_prev;
    iter = 0;
    do {
        x_prev = x;
        double fx, dfx;
        if (use_round) {
            fx = Round(f(x_prev), Delta);
            dfx = Round(df(x_prev), Delta);
        } else {
            fx = f(x_prev);
            dfx = df(x_prev);
        }
        if (fabs(dfx) < 1e-12) {
            cerr << "Ошибка: производная близка к нулю" << endl;
            break;
        }
        x = x_prev - fx / dfx;
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
    
    cout << "=== Исследование 1: зависимость числа итераций от eps ===" << endl;
    cout << "eps\t\tкорень\t\tитераций" << endl;
    cout << fixed << setprecision(12); 
    for (double eps = 0.1; eps >= 1e-6; eps /= 10) {
        double root = newton(x0, eps, iter, false);
        cout << eps << "\t\t" << root << "\t\t" << iter << endl;
    }
    
    cout << "\n=== Исследование 2: влияние округления (eps = 1e-6) ===" << endl;
    cout << "Delta\t\tкорень\t\tитераций" << endl;
    double eps_fixed = 1e-6;
    double Delta_values[] = {0.0, 0.001, 0.01, 0.1, 1.0};
    for (int i = 0; i < 5; i++) {
        Delta = Delta_values[i];
        double root = newton(x0, eps_fixed, iter, true);
        cout << Delta << "\t\t" << root << "\t\t" << iter << endl;
    }
    
    return 0;
}
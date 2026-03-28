#include <iostream>
#include <iomanip>
#include <cmath>
#include <clocale>

using namespace std;

double DELTA = 0.0;

double Round(double, double);
double ITER(double, double, int&);

int f_calls = 0;

// Функция φ(x)
double F(double x) {
    f_calls++;
    if (f_calls > 1000) return x; // Защита от зацикливания при грубых округлениях
    
    double y = exp(1.0 / sqrt(x));
    if (DELTA > 1e-9) return Round(y, DELTA);
    return y;
}

// Производная φ'(x) (нужна для отчета, в самом ITER не вызывается)
double F1(double x) {
    double y = -exp(1.0 / sqrt(x)) / (2.0 * x * sqrt(x));
    if (DELTA > 1e-9) return Round(y, DELTA);
    return y;
}

int main() {
    setlocale(LC_ALL, "Russian");
    
    double x0 = 2.0; 
    int N;

    cout << "Лаб. 4: Метод простых итераций\n";
    cout << "Уравнение: ln(x)^2 - 1/x = 0  =>  x = exp(1/sqrt(x))\n";
    cout << "Начальное приближение x0: " << x0 << "\n\n";

    DELTA = 0.0;
    
    cout << "6.1) Зависимость итераций (N) от точности (Eps):\n";
    cout << fixed << setprecision(10);
    cout << setw(10) << "Eps" << setw(18) << "X" << setw(10) << "N\n";
    
    double exactX = 0;
    for (double eps = 0.1; eps >= 0.9e-6; eps /= 10) {
        f_calls = 0;
        double x = ITER(x0, eps, N);
        if (f_calls > 1000) N = -1; // Маркер зацикливания
        cout << setw(10) << setprecision(6) << eps << setw(18) << setprecision(9) << x << setw(10) << N << endl;
        if (eps < 2e-6) exactX = x;
    }

    cout << "\n6.2) Чувствительность к ошибкам (Eps=1e-6):\n";
    cout << setw(10) << "Delta" << setw(18) << "X" << setw(10) << "N" << setw(18) << "Error\n";
    
    double deltas[] = {0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 0.0};
    
    for (double d : deltas) {
        DELTA = d;
        f_calls = 0;
        double x = ITER(x0, 1e-6, N);
        if (f_calls > 1000) {
            cout << setw(10) << setprecision(6) << d << setw(18) << "Зацикливание" << setw(10) << ">1000" << setw(18) << "---" << endl;
        } else {
            cout << setw(10) << setprecision(6) << d << setw(18) << setprecision(9) << x << setw(10) << N << setw(18) << setprecision(9) << abs(x - exactX) << endl;
        }
    }

    return 0;
}

#include "methods.cpp"

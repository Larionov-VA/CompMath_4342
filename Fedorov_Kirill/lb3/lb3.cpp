#define __NEWTON
#include <iostream>
#include <iomanip>
#include <cmath>
#include <clocale>

using namespace std;

// Глобальная переменная для погрешности (пункт 5)
double DELTA = 0.0;

// Объявления функций из methods.cpp
double Round(double, double);
double NEWTON(double, double, int&);

// Наша функция: ln(x)^2 - 1/x
double F(double x) {
    double y = pow(log(x), 2) - 1.0 / x;
    if (DELTA > 1e-9) return Round(y, DELTA);
    return y;
}

// Первая производная: 2*ln(x)/x + 1/x^2 = (2*x*ln(x) + 1)/x^2
double F1(double x) {
    double y = (2.0 * x * log(x) + 1.0) / (x * x);
    if (DELTA > 1e-9) return Round(y, DELTA);
    return y;
}

int main() {
    setlocale(LC_ALL, "Russian");
    
    double a = 2.0, b = 2.5; 
    double x0 = 2.0; // Т.к. f(2)*f''(2) > 0
    int N;

    cout << "Лаб. 3: Метод Ньютона (касательных)\n";
    cout << "Уравнение: ln(x)^2 - 1/x = 0\n";
    cout << "Отрезок: [" << a << ", " << b << "]\n";
    cout << "Начальное приближение x0: " << x0 << "\n\n";

    DELTA = 0.0;
    
    // Пункт 5.1: Зависимость от Eps
    cout << "5.1) Зависимость итераций (N) от точности (Eps):\n";
    cout << fixed << setprecision(10);
    cout << setw(10) << "Eps" << setw(18) << "X" << setw(10) << "N\n";
    
    double exactX = 0;
    for (double eps = 0.1; eps >= 0.9e-6; eps /= 10) {
        double x = NEWTON(x0, eps, N);
        cout << setw(10) << setprecision(6) << eps << setw(18) << setprecision(9) << x << setw(10) << N << endl;
        if (eps < 2e-6) exactX = x;
    }

    // Пункт 5.2: Чувствительность к ошибкам
    cout << "\n5.2) Чувствительность к ошибкам (Eps=1e-6):\n";
    cout << setw(10) << "Delta" << setw(18) << "X" << setw(10) << "N" << setw(18) << "Error\n";
    
    double deltas[] = {0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 0.0};
    
    for (double d : deltas) {
        DELTA = d;
        // Перехват исключений не реализован в NEWTON, 
        // поэтому при грубых округлениях производная может стать 0 и программа завершится.
        // Но попробуем выполнить как есть:
        double x = NEWTON(x0, 1e-6, N);
        cout << setw(10) << setprecision(6) << d << setw(18) << setprecision(9) << x << setw(10) << N << setw(18) << setprecision(9) << abs(x - exactX) << endl;
    }

    return 0;
}

#include "methods.cpp"

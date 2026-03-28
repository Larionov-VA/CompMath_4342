#include <iostream>
#include <iomanip>
#include <cmath>
#include <clocale>

using namespace std;

double DELTA = 0.0;

double Round(double, double);
double HORDA(double, double, double, int&);

// Наша функция: ln(x)^2 - 1/x
double F(double x) {
    double y = pow(log(x), 2) - 1.0 / x;
     // Если DELTA задана (>0), то с ошибкой, иначе точно
    if (DELTA > 1e-9) return Round(y, DELTA);
    return y;
}

int main() {
    setlocale(LC_ALL, "Russian");
    
    double a = 2.0, b = 2.5; 
    int N;

    cout << "Лаб. 2: Метод хорд\n";
    cout << "Уравнение: ln(x)^2 - 1/x = 0\n";
    cout << "Отрезок: [" << a << ", " << b << "]\n\n";

    // Проверка корректности отрезка
    DELTA = 0.0;
    if (F(a) * F(b) >= 0) {
        cout << "Ошибка: на концах отрезка функция одного знака!\n";
        return 1;
    }
    
    // Пункт 4: Зависимость от Eps
    cout << "4) Зависимость итераций (N) от точности (Eps):\n";
    cout << fixed << setprecision(10);
    cout << setw(10) << "Eps" << setw(18) << "X" << setw(10) << "N\n";
    
    double exactX = 0;
    for (double eps = 0.1; eps >= 0.9e-6; eps /= 10) {
        double x = HORDA(a, b, eps, N);
        cout << setw(10) << setprecision(6) << eps << setw(18) << setprecision(9) << x << setw(10) << N << endl;
        if (eps < 2e-6) exactX = x;
    }

    // Пункт 5: Чувствительность к ошибкам
    cout << "\n5) Чувствительность к ошибкам (Eps=1e-6):\n";
    cout << setw(10) << "Delta" << setw(18) << "X" << setw(10) << "N" << setw(18) << "Error\n";
    
    double deltas[] = {0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 0.0};
    
    for (double d : deltas) {
        DELTA = d;
        double x = HORDA(a, b, 1e-6, N);
        cout << setw(10) << setprecision(6) << d << setw(18) << setprecision(9) << x << setw(10) << N << setw(18) << setprecision(9) << abs(x - exactX) << endl;
    }

    return 0;
}

#include "methods.cpp"

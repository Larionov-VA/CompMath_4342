#include <iostream>
#include <cmath>
#include <iomanip>
#include "methods.h"

using namespace std;   

bool use_rounding = false;
double current_delta = 0.1; 

double F(double x) {
    double value = atan(x) - 1.0/x;

    if (use_rounding) {
        value = Round(value, current_delta);
    }
    return value;
}

int main(){
    double a = 1.0;
    double b = 2.0;
    int N = 0;

    use_rounding = false; 
    cout << "--- Метод Хорд: Зависимость итераций от Eps ---" << endl;
    cout << "   Eps      |   Корень x   | Итераций N" << endl;
    cout << "------------------------------------------" << endl;

    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        double root = HORDA(a, b, eps, N);
        cout << fixed << setprecision(6) << eps << "    | " << root << "     | " << N << endl;
    }

    use_rounding = true; 
    cout << "\nТаблица 2 – Чувствительность метода хорд к ошибкам (Delta)" << endl;
    cout << "| Eps      | Delta    | Результат метода |" << endl;
    cout << "|----------|----------|------------------|" << endl;

    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        for (double delta = 0.1; delta >= 0.0000009; delta /= 10.0) {
            current_delta = delta;
            double root = HORDA(a, b, eps, N);
            cout << "| " << fixed << setprecision(6) << eps << " | " 
                 << delta << " | " << root << " |" << endl;
        }
        cout << "|----------|----------|------------------|" << endl;
    }

    return 0;
}
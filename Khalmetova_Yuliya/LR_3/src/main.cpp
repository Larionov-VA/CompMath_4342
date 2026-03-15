#define __NEWTON
#include <iostream>
#include <cmath>
#include <iomanip>
#include "methods.h"

using namespace std;

bool use_rounding = false;
double current_delta = 0.1;

double F(double x) {
    double value = atan(x) - 1.0 / x;
    if (use_rounding) {
        value = Round(value, current_delta);
    }
    return value;
}

double F1(double x) {
    double value = (1.0 / (1.0 + x * x)) + (1.0 / (x * x));
    if (use_rounding) {
        value = Round(value, current_delta);
    }
    return value;
}

int main() {
    double x0 = 1.0; 
    int N = 0;
    
    cout << "--- Часть 1: Зависимость количества итераций от точности Eps ---" << endl;
    cout << "   Eps   | Модиф. Eps0 |   Корень x   | Итераций N" << endl;
    cout << "------------------------------------------------------" << endl;
    
    use_rounding = false; 
    
    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        double eps0 = 0.6 * sqrt(eps);
        double root = NEWTON(x0, eps0, N);
        
        cout << fixed << setprecision(6) << eps << " | " 
             << eps0 << "    | " << root << "     | " << N << endl;
    }

    cout << "\n--- Часть 2: Исследование чувствительности метода к ошибкам (Delta) ---" << endl;
    cout << "| Eps          | Модиф. Eps0 |  Delta   | Результат метода |" << endl;
    cout << "|--------------|-------------|----------|------------------|" << endl;
    
    use_rounding = true;
    
    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        double eps0 = 0.6 * sqrt(eps); 
        
        for (double delta = 0.1; delta >= 0.0000009; delta /= 10.0) {
            current_delta = delta;
            
            double root = NEWTON(x0, eps0, N);
            
            cout << "| " << fixed << setprecision(6) << eps << "     | "
                 << eps0 << "    | " << delta << " | " << root << " |" << endl;
        }
        cout << "|--------------|-------------|----------|------------------|" << endl;
    }

    return 0;
}
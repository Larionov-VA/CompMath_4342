#include <iostream>
#include <cmath>
#include <iomanip>
#include "methods.h" 

using namespace std;

bool use_rounding = false;
double current_delta = 0.1;
int iteration_guard = 0;

double F(double x) {
    iteration_guard++;
    
    if (iteration_guard > 1000) {
        return x;
    }

    double value = 1.0 / atan(x); 
    
    if (use_rounding) {
        value = Round(value, current_delta);
    }
    return value;
}

int main() {
    double x0 = 1.0; 
    int N = 0;

    use_rounding = false;
    cout << "--- Метод простых итераций: Зависимость итераций от Eps ---" << endl;
    cout << "   Eps      |   Корень x   | Итераций N" << endl;
    cout << "------------------------------------------" << endl;
    
    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        iteration_guard = 0;
        double root = ITER(x0, eps, N);
        cout << fixed << setprecision(6) << eps << "  | " << root << "     | " << N << endl;
    }

    use_rounding = true;
    cout << "\nТаблица 2 - Чувствительность метода итераций к ошибкам (Delta)" << endl;
    cout << "| Eps      | Delta    | Результат метода | Итераций (N) |" << endl;
    cout << "|----------|----------|------------------|--------------|" << endl;
    
    for (double eps = 0.1; eps >= 0.0000009; eps /= 10.0) {
        for (double delta = 0.1; delta >= 0.0000009; delta /= 10.0) {
            current_delta = delta;
            iteration_guard = 0;
            
            double root = ITER(x0, eps, N);
            
            if (iteration_guard > 1000) {
                cout << "| " << fixed << setprecision(6) << eps << " | " 
                     << delta << " | " << root << " |  Зациклился  |" << endl;
            } else {
                cout << "| " << fixed << setprecision(6) << eps << " | " 
                     << delta << " | " << root << " | " << setw(12) << N << " |" << endl;
            }
        }
        cout << "|----------|----------|------------------|--------------|" << endl;
    }

    return 0;
}
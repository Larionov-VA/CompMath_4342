#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

// 1. Подынтегральная функция для варианта №17
// f(x) = (x / (1 + x^2)) * exp(-x^2)
double f(double x) {
    return (x / (1.0 + x * x)) * exp(-x * x);
}

// 2. Составная формула прямоугольников (средних)
double rectangleRule(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double x_mid = a + i * h + h / 2.0;
        sum += f(x_mid);
    }
    return h * sum;
}

// 3. Составная формула трапеций
double trapezoidRule(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }
    return h * sum;
}

// 4. Составная формула Симпсона
double simpsonRule(double a, double b, int n) {
    if (n % 2 != 0) n++; 
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x_i = a + i * h;
        if (i % 2 != 0)
            sum += 4.0 * f(x_i);
        else
            sum += 2.0 * f(x_i);
    }
    return (h / 3.0) * sum;
}

void solveForEps(double a, double b, double eps) {
    cout << "\n>>> Точность eps = " << eps << " <<<" << endl;
    cout << left << setw(18) << "Метод" << " | " << setw(6) << "n" 
         << " | " << setw(14) << "Значение Ih2" << " | " << setw(14) << "Уточненное I" << " | " << "Погр. (Рунге)" << endl;
    cout << string(85, '-') << endl;

    string names[] = {"Прямоугольников", "Трапеций", "Симпсона"};
    int divisors[] = {3, 3, 15};

    for (int k = 0; k < 3; ++k) {
        int n = 2;
        double Ih, Ih2, error, final_I;
        while (true) {
            if (k == 0) { Ih = rectangleRule(a, b, n); Ih2 = rectangleRule(a, b, 2 * n); }
            else if (k == 1) { Ih = trapezoidRule(a, b, n); Ih2 = trapezoidRule(a, b, 2 * n); }
            else { Ih = simpsonRule(a, b, n); Ih2 = simpsonRule(a, b, 2 * n); }

            error = abs(Ih2 - Ih) / (double)divisors[k];
            if (error < eps) {
                double diff = abs(Ih2 - Ih);
                if (k == 0) final_I = Ih2 + diff / 3.0;
                else if (k == 1) final_I = Ih2 - diff / 3.0;
                else final_I = Ih2 - diff / 15.0;
                
                cout << left << setw(18) << names[k] << " | " << setw(6) << 2 * n 
                     << " | " << fixed << setprecision(8) << setw(14) << Ih2 
                     << " | " << setw(14) << final_I 
                     << " | " << scientific << setprecision(2) << error << endl;
                break;
            }
            n *= 2;
            if (n > 2000000) break;
        }
    }
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double precisions[] = {0.01, 0.001, 0.0001};

    cout << "Лабораторная работа N6. Вариант 17." << endl;
    cout << "Интеграл x/(1+x^2)*exp(-x^2) от 0 до 1" << endl;

    for (double eps : precisions) {
        solveForEps(a, b, eps);
    }
    return 0;
}
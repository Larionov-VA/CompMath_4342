#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

// Константа Пи
const double PI = acos(-1.0);

// 1) Программа-функция для вычисления подынтегральной функции.
// f(x) = sqrt(x) * exp(-x^2)
double f(double x) {
    return sqrt(x) * exp(-x * x);
}

// Формула средних прямоугольников
double rectangleRule(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        double x_mid = a + i * h + h / 2.0;
        sum += f(x_mid);
    }

    return h * sum;
}

// Формула трапеций
double trapezoidalRule(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = (f(a) + f(b)) / 2.0;

    for (int i = 1; i < n; ++i) {
        sum += f(a + i * h);
    }

    return h * sum;
}

// Формула Симпсона
double simpsonRule(double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; ++i) {
        int factor = (i % 2 == 0) ? 2 : 4;
        sum += factor * f(a + i * h);
    }

    return h * sum / 3.0;
}

// Оценка по Рунге и поиск нужного n
void calculateWithRunge(const string& methodName,
                        double (*method)(double, double, int),
                        double a, double b, double eps, int rungeDivisor)
{
    int n = 2;
    double I_h = method(a, b, n);
    double I_h2;
    double error;
    double rungeTerm;

    do {
        n *= 2;
        I_h2 = method(a, b, n);

        rungeTerm = (I_h2 - I_h) / rungeDivisor;
        error = abs(rungeTerm);

        I_h = I_h2;
    } while (error > eps);

    double refined_I = I_h2 + rungeTerm;

    cout << left << setw(18) << methodName
         << " | n = " << setw(5) << n
         << " | I_h/2 = " << fixed << setprecision(8) << I_h2
         << " | Уточн. I = " << refined_I
         << " | Погрешность: " << scientific << setprecision(2) << error
         << fixed << endl;
}

int main() {
    // Границы интегрирования
    double a = PI / 2.0;
    double b = PI;

    int count;
    cout << "Введите количество значений epsilon: ";
    cin >> count;

    vector<double> eps_values(count);

    cout << "Введите " << count << " значений epsilon через пробел: ";
    for (int i = 0; i < count; ++i) {
        cin >> eps_values[i];
    }

    cout << "\nИнтегрирование функции f(x) = sqrt(x) * exp(-x^2)\n";
    cout << "На отрезке [pi/2, pi]\n";
    cout << string(95, '=') << endl;

    // Цикл по всем введенным значениям epsilon
    for (double eps : eps_values) {
        cout << "ТРЕБУЕМАЯ ТОЧНОСТЬ (Epsilon) = " << defaultfloat << eps << endl;
        cout << string(95, '-') << endl;

        // Для прямоугольников и трапеций p = 2, поэтому 2^2 - 1 = 3
        calculateWithRunge("Прямоугольников", rectangleRule, a, b, eps, 3);
        calculateWithRunge("Трапеций", trapezoidalRule, a, b, eps, 3);

        // Для Симпсона p = 4, поэтому 2^4 - 1 = 15
        calculateWithRunge("Симпсона", simpsonRule, a, b, eps, 15);

        cout << string(95, '=') << endl;
    }

    return 0;
}
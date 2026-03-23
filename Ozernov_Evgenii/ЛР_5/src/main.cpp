#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include "methods.h"

// f(x) = arcsin(2x/(1+x^2)) - e^(-x^2)
double f(double x) {
    double t = (2.0 * x) / (1.0 + x * x);

    // защита от погрешностей double: зажимаем в [-1;1]
    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;

    return std::asin(t) - std::exp(-x * x);
}

// f'(x) = t'(x)/sqrt(1-t(x)^2) + 2x*e^(-x^2)
// t(x) = 2x/(1+x^2)
// t'(x) = 2(1 - x^2)/(1 + x^2)^2
double fp(double x) {
    double denom = (1.0 + x * x);
    double t = (2.0 * x) / denom;

    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;

    double tp = 2.0 * (1.0 - x * x) / (denom * denom);
    double s = std::sqrt(std::max(0.0, 1.0 - t * t));

    // если s очень мал (рядом с t=±1), защищаемся
    double part1 = (s < 1e-15) ? 0.0 : (tp / s);
    double part2 = 2.0 * x * std::exp(-x * x);

    return part1 + part2;
}

// Численная f''(x) для выбора x0 по условию f(x0)*f''(x0)>0
double f2_num(double x) {
    double h = 1e-5 * std::max(1.0, std::abs(x));
    return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h);
}

// Авто-отделение корня: ищем смену знака на сетке
bool find_bracket(double x_min, double x_max, double step, double &Left, double &Right) {
    double prev = x_min;
    double fprev = f(prev);

    for (double x = x_min + step; x <= x_max; x += step) {
        double fx = f(x);

        if (fprev == 0.0) { Left = prev; Right = prev; return true; }
        if (fprev * fx < 0.0) { Left = prev; Right = x; return true; }

        prev = x;
        fprev = fx;
    }
    return false;
}

int main() {
    std::cout << std::fixed << std::setprecision(12);

    // 1) Отделяем корень (это только диапазон поиска!)
    double Left = 0.0, Right = 0.0;
    if (!find_bracket(-5.0, 5.0, 0.1, Left, Right)) {
        std::cerr << "Не удалось отделить корень на [-5;5] с шагом 0.1\n";
        return 1;
    }

    std::cout << "Отделение корня (авто): [" << Left << "; " << Right << "]\n";
    std::cout << "f(Left)=" << f(Left) << ", f(Right)=" << f(Right) << "\n";

    // 2) Выбор x0: требуется f(x0)*f''(x0) > 0
    double x0 = 0.5 * (Left + Right);
    double condL = f(Left)  * f2_num(Left);
    double condR = f(Right) * f2_num(Right);
    double condM = f(x0)    * f2_num(x0);

    if (condL > 0.0) x0 = Left;
    else if (condR > 0.0) x0 = Right;
    else if (condM > 0.0) x0 = 0.5 * (Left + Right);
    // если ни одна точка не подошла (редко), оставляем середину

    std::cout << "Выбрано начальное приближение x0=" << x0
              << " (проверка f(x0)*f''(x0)>0)\n\n";

    // 3) Опорное значение корня (эталон): Delta=0, очень высокая точность
    auto ref = NEWTON(f, fp, Left, Right, x0, 1e-14, 1000000, 0.0);
    if (!ref.ok) {
        std::cerr << "Не удалось получить reference-корень методом Ньютона\n";
        return 2;
    }
    std::cout << "Reference root (Eps=1e-14, Delta=0): x=" << ref.x
              << "  iters=" << ref.iters << "\n\n";

    // 4) Скорость сходимости: Iterations vs Eps (0.1 .. 1e-6)
    std::ofstream feps("iterations_vs_eps_newton.csv");
    feps << "Eps;Iterations;Root;LastStep\n";

    for (int i = 0; i <= 50; ++i) {
        double t = i / 50.0;          // 0..1
        double p = -1.0 - 5.0 * t;    // -1..-6
        double Eps = std::pow(10.0, p);

        auto r = NEWTON(f, fp, Left, Right, x0, Eps, 1000000, 0.0);

        feps << std::setprecision(16)
             << Eps << ";" << r.iters << ";" << r.x << ";" << r.last_step << "\n";
    }
    feps.close();
    std::cout << "CSV сохранён: iterations_vs_eps_newton.csv\n";

    // 5) Обусловленность: влияние Delta (округление f и f')
    const double EpsFix = 1e-6;
    std::vector<double> deltas = {0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};

    std::ofstream fdel("sensitivity_vs_delta_newton.csv");
    fdel << "Delta;Eps;Iterations;Root;AbsErrorToRef\n";

    for (double Delta : deltas) {
        auto r = NEWTON(f, fp, Left, Right, x0, EpsFix, 1000000, Delta);
        double err = std::abs(r.x - ref.x);

        fdel << std::setprecision(16)
             << Delta << ";" << EpsFix << ";" << r.iters << ";" << r.x << ";" << err << "\n";
    }
    fdel.close();
    std::cout << "CSV сохранён: sensitivity_vs_delta_newton.csv\n";

    return 0;
}
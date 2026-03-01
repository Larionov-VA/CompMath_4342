#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include "methods.h"

// f(x) = arcsin(2x/(1+x^2)) - e^(-x^2)
double f(double x) {
    double t = (2.0 * x) / (1.0 + x * x);

    // защита от погрешностей double
    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;

    return std::asin(t) - std::exp(-x * x);
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

    // 1) Отделение корня
    double Left = 0.0, Right = 0.0;
    if (!find_bracket(-5.0, 5.0, 0.1, Left, Right)) {
        std::cerr << "Не удалось отделить корень на [-5,5] с шагом 0.1\n";
        return 1;
    }

    std::cout << "Отделение корня (авто): [" << Left << "; " << Right << "]\n";
    std::cout << "f(Left)=" << f(Left) << ", f(Right)=" << f(Right) << "\n\n";

    // 2) Опорный (reference) корень: высокая точность, без округления
    auto ref = HORDA(f, Left, Right, 1e-14, 2000000, 0.0);
    if (!ref.ok) {
        std::cerr << "Не удалось получить reference-корень методом хорд\n";
        return 2;
    }
    std::cout << "Reference root (Eps=1e-14, Delta=0): x=" << ref.x
              << "  iters=" << ref.iters << "\n\n";

    // 3) Скорость сходимости: Iterations vs Eps (0.1 .. 1e-6)
    std::ofstream feps("iterations_vs_eps_horda.csv");
    feps << "Eps;Iterations;Root;IntervalLen\n";

    for (int i = 0; i <= 50; ++i) {
        double t = i / 50.0;          // 0..1
        double p = -1.0 - 5.0 * t;    // -1..-6
        double Eps = std::pow(10.0, p);

        auto r = HORDA(f, Left, Right, Eps, 2000000, 0.0);
        double len = std::abs(r.b - r.a);

        feps << std::setprecision(16)
             << Eps << ";" << r.iters << ";" << r.x << ";" << len << "\n";
    }
    feps.close();
    std::cout << "CSV сохранён: iterations_vs_eps_horda.csv\n";

    // 4) Обусловленность: влияние Delta (округление значений функции)
    const double EpsFix = 1e-6;
    std::vector<double> deltas = {0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};

    std::ofstream fdel("sensitivity_vs_delta_horda.csv");
    fdel << "Delta;Eps;Iterations;Root;AbsErrorToRef\n";

    for (double Delta : deltas) {
        auto r = HORDA(f, Left, Right, EpsFix, 2000000, Delta);
        double err = std::abs(r.x - ref.x);

        fdel << std::setprecision(16)
             << Delta << ";" << EpsFix << ";" << r.iters << ";" << r.x << ";" << err << "\n";
    }
    fdel.close();
    std::cout << "CSV сохранён: sensitivity_vs_delta_horda.csv\n";

    return 0;
}
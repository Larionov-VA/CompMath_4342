#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "methods.h"

// ============================================================
// ВАРИАНТ 18: f(x) = e^(1/x^2) - ln(x)
// Область определения: x > 0
// Интервал изоляции: [3; 4]
// ============================================================
double f(double x) {
    if (x <= 0.0) {
        throw std::runtime_error("Function is undefined for x <= 0");
    }
    return std::exp(1.0 / (x * x)) - std::log(x);
}

// ============================================================
// Головная программа
// ============================================================
int main() {
    const double Left = 3.0;
    const double Right = 4.0;
    const int Nmax = 2000000;

    std::cout << std::fixed << std::setprecision(12);
    std::cout << "Laboratornaya rabota #4. Metod hord.\n";
    std::cout << "Variant 18: f(x) = exp(1/x^2) - ln(x)\n";
    std::cout << "Interval [" << Left << "; " << Right << "]\n";
    std::cout << "f(Left)  = " << f(Left) << "\n";
    std::cout << "f(Right) = " << f(Right) << "\n\n";

    const auto reference = HORDA(f, Left, Right, 1e-14, Nmax, 0.0);
    if (!reference.ok) {
        std::cerr << "Reference root was not obtained.\n";
        return 1;
    }

    std::cout << "Reference root (eps = 1e-14, delta = 0):\n";
    std::cout << "x* = " << reference.x
              << ", iterations = " << reference.iters << "\n\n";

    std::vector<double> eps_values;
    for (int i = 0; i <= 40; ++i) {
        const double alpha = static_cast<double>(i) / 40.0;
        // Изменяем eps от 1e-1 до 1e-6 по логарифмической шкале.
        const double power = -1.0 - 5.0 * alpha;
        eps_values.push_back(std::pow(10.0, power));
    }

    // ============================================================
    // Исследование зависимости числа итераций от точности eps
    // ============================================================
    std::ofstream eps_file("eps_iter_profile.csv");
    eps_file << "Eps;Iterations;Root;IntervalLen;Ok\n";

    std::cout << "Speed of convergence study:\n";
    for (double eps : eps_values) {
        const auto result = HORDA(f, Left, Right, eps, Nmax, 0.0);
        const double interval_len = std::abs(result.b - result.a);

        eps_file << std::setprecision(16)
                 << eps << ';'
                 << result.iters << ';'
                 << result.x << ';'
                 << interval_len << ';'
                 << (result.ok ? 1 : 0) << '\n';

        std::cout << "eps = " << std::setw(12) << eps
                  << "  iterations = " << std::setw(3) << result.iters
                  << "  x = " << result.x << '\n';
    }
    eps_file.close();

    // ============================================================
    // Исследование чувствительности к ошибкам округления delta
    // ============================================================
    const double fixed_eps = 1e-8;
    const std::vector<double> delta_values = {
        0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8
    };

    std::ofstream delta_file("delta_impact_profile.csv");
    delta_file << "Delta;Eps;Iterations;Root;AbsErrorToRef;Ok\n";

    std::cout << "\nSensitivity study (fixed eps = 1e-8):\n";
    for (double delta : delta_values) {
        const auto result = HORDA(f, Left, Right, fixed_eps, Nmax, delta);
        const double error = result.ok
            ? std::abs(result.x - reference.x)
            : std::numeric_limits<double>::quiet_NaN();

        delta_file << std::setprecision(16)
                   << delta << ';'
                   << fixed_eps << ';'
                   << result.iters << ';'
                   << result.x << ';'
                   << error << ';'
                   << (result.ok ? 1 : 0) << '\n';

        std::cout << "delta = " << std::setw(10) << delta
                  << "  ok = " << result.ok
                  << "  iterations = " << std::setw(3) << result.iters
                  << "  x = " << result.x << '\n';
    }
    delta_file.close();

    return 0;
}

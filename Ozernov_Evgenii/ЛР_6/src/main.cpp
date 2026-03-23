#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include "methods.h"

// f(x) = arcsin(2x/(1+x^2)) - e^(-x^2)
double f(double x) {
    double t = (2.0 * x) / (1.0 + x * x);

    // защита от погрешностей double: зажимаем в [-1;1]
    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;

    return std::asin(t) - std::exp(-x * x);
}

// Выбранное преобразование:
// 2*atan(x) = e^(-x^2)  =>  x = tan(e^(-x^2)/2)
double phi(double x) {
    return std::tan(0.5 * std::exp(-x * x));
}

// phi'(x) = -x * e^(-x^2) * sec^2(e^(-x^2)/2)
double dphi(double x) {
    double a = 0.5 * std::exp(-x * x);
    double c = std::cos(a);
    double sec2 = 1.0 / (c * c);
    return -x * std::exp(-x * x) * sec2;
}

// Авто-отделение корня: ищем смену знака на сетке
bool find_bracket(double x_min, double x_max, double step, double &Left, double &Right) {
    double prev = x_min;
    double fprev = f(prev);

    for (double x = x_min + step; x <= x_max + 1e-15; x += step) {
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

    // 1) Отделяем корень
    double Left = 0.0, Right = 0.0;
    if (!find_bracket(-5.0, 5.0, 0.1, Left, Right)) {
        std::cerr << "Не удалось отделить корень на [-5;5] с шагом 0.1\n";
        return 1;
    }

    std::cout << "Отделение корня (авто): [" << Left << "; " << Right << "]\n";
    std::cout << "f(Left)=" << f(Left) << ", f(Right)=" << f(Right) << "\n\n";

    // 2) Проверяем отображение интервала в себя и условие сходимости
    double phiL = phi(Left);
    double phiR = phi(Right);
    double q = estimate_q(dphi, Left, Right, 5000, 0.0);

    std::cout << "Выбрана функция phi(x) = tan(exp(-x^2)/2)\n";
    std::cout << "phi(Left)=" << phiL << ", phi(Right)=" << phiR << "\n";
    std::cout << "q = max |phi'(x)| на [Left;Right] = " << q << "\n";
    std::cout << "Так как q < 1, метод простых итераций сходится.\n\n";

    // 3) Выбор x0: любая точка из [Left;Right], берем середину
    double x0 = 0.5 * (Left + Right);
    std::cout << "Выбрано начальное приближение x0=" << x0 << "\n\n";

    // 4) Опорное значение корня (эталон): Delta=0, очень высокая точность
    auto ref = ITER(phi, dphi, Left, Right, x0, 1e-14, 1000000, 0.0);
    if (!ref.ok) {
        std::cerr << "Не удалось получить reference-корень методом простых итераций\n";
        return 2;
    }
    std::cout << "Reference root (Eps=1e-14, Delta=0): x=" << ref.x
              << "  iters=" << ref.iters << "\n";
    std::cout << "Проверка f(x_ref)=" << f(ref.x) << "\n\n";

    // 5) Скорость сходимости: Iterations vs Eps (0.1 .. 1e-10)
    std::ofstream feps("iterations_vs_eps_iter.csv");
    feps << "Eps;Iterations;Root;LastStep\n";

    for (int i = 0; i <= 90; ++i) {
        double t = i / 90.0;          // 0..1
        double p = -1.0 - 9.0 * t;    // -1..-10
        double Eps = std::pow(10.0, p);

        auto r = ITER(phi, dphi, Left, Right, x0, Eps, 1000000, 0.0);

        feps << std::setprecision(16)
             << Eps << ";" << r.iters << ";" << r.x << ";" << r.last_step << "\n";
    }
    feps.close();
    std::cout << "CSV сохранён: iterations_vs_eps_iter.csv\n";

    // 6) Обусловленность: влияние Delta (округление phi и phi')
    const double EpsFix = 1e-6;
    std::vector<double> deltas = {0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};

    std::ofstream fdel("sensitivity_vs_delta_iter.csv");
    fdel << "Delta;Eps;Iterations;Root;AbsErrorToRef\n";

    for (double Delta : deltas) {
        auto r = ITER(phi, dphi, Left, Right, x0, EpsFix, 1000000, Delta);
        double err = std::abs(r.x - ref.x);

        fdel << std::setprecision(16)
             << Delta << ";" << EpsFix << ";" << r.iters << ";" << r.x << ";" << err << "\n";
    }
    fdel.close();
    std::cout << "CSV сохранён: sensitivity_vs_delta_iter.csv\n\n";

    // 7) Один рабочий запуск для вывода результата в консоль
    const double EpsWork = 1e-8;
    const double DeltaWork = 1e-10;
    auto work = ITER(phi, dphi, Left, Right, x0, EpsWork, 1000000, DeltaWork);

    std::cout << "Рабочий запуск: Eps=" << EpsWork << ", Delta=" << DeltaWork << "\n";
    std::cout << "Корень x=" << work.x << "\n";
    std::cout << "f(x)=" << f(work.x) << "\n";
    std::cout << "Итераций=" << work.iters << "\n";
    std::cout << "Последний шаг=" << work.last_step << "\n";
    std::cout << "Абсолютная ошибка относительно эталона=" << std::abs(work.x - ref.x) << "\n";

    // Так как phi'(x) < 0 на [Left;Right], по методичке погрешность не превышает Eps
    double sign_check = dphi(0.5 * (Left + Right));
    if (sign_check < 0.0) {
        std::cout << "Так как phi'(x) < 0 на интервале, теоретическая оценка погрешности: <= Eps\n";
    } else {
        std::cout << "Теоретическая оценка погрешности: q*Eps/(1-q) = "
                  << q * EpsWork / (1.0 - q) << "\n";
    }

    return 0;
}

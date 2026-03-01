#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include "methods.h"

// f(x) = arcsin(2x/(1+x^2)) - e^(-x^2)
double f(double x) {
    double t = (2.0 * x) / (1.0 + x * x);

    // на всякий случай зажмём в [-1,1] (защита от погрешностей double)
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

    // 1) Отделяем корень
    double Left = 0.0, Right = 0.0;
    if (!find_bracket(-5.0, 5.0, 0.1, Left, Right)) {
        std::cerr << "Не удалось отделить корень на [-5,5] с шагом 0.1\n";
        return 1;
    }

    std::cout << "Отделение корня (авто): [" << Left << ", " << Right << "]\n";
    std::cout << "f(Left)=" << f(Left) << ", f(Right)=" << f(Right) << "\n\n";

    // 2) Находим опорное "почти точное" значение корня (для сравнения)
    auto ref = BISECT(f, Left, Right, 1e-14, 2000000, 0.0);
    if (!ref.ok) {
        std::cerr << "Не удалось получить reference-корень\n";
        return 2;
    }
    std::cout << "Reference root (Eps=1e-14, Delta=0): x=" << ref.x
              << "  iters=" << ref.iters << "\n\n";

    // 3) Зависимость итераций от Eps: от 1e-1 до 1e-6 (лог-шкала, 51 точка)
    std::ofstream feps("iterations_vs_eps.csv");
    feps << "Eps,Iterations,Root,IntervalLen\n";

    std::cout << "Таблица (часть) Iterations vs Eps:\n";
    std::cout << "Eps\t\titers\troot\n";

    for (int i = 0; i <= 50; ++i) {
        double t = i / 50.0;             // 0..1
        double p = -1.0 - 5.0 * t;       // -1..-6
        double Eps = std::pow(10.0, p);

        auto r = BISECT(f, Left, Right, Eps, 1000000, 0.0);
        double len = std::abs(r.b - r.a);

        feps << std::setprecision(16) << Eps << ","
             << r.iters << ","
             << r.x << ","
             << len << "\n";

        if (i % 10 == 0 || i == 50) {
            std::cout << std::setprecision(6) << Eps << "\t"
                      << std::setprecision(0) << r.iters << "\t"
                      << std::setprecision(12) << r.x << "\n";
        }
    }
    feps.close();
    std::cout << "\nCSV сохранён: iterations_vs_eps.csv\n\n";

    // 4) Чувствительность к ошибкам (округление f(...) с точностью Delta)
    //    фиксируем Eps=1e-6, меняем Delta
    const double EpsFix = 1e-6;
    std::vector<double> deltas = {0.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6};

    std::ofstream fdel("sensitivity_vs_delta.csv");
    fdel << "Delta,Eps,Iterations,Root,AbsErrorToRef\n";

    std::cout << "Чувствительность (Eps=1e-6):\n";
    std::cout << "Delta\t\titers\troot\t\t\t|x-ref|\n";

    for (double Delta : deltas) {
        auto r = BISECT(f, Left, Right, EpsFix, 1000000, Delta);
        double err = std::abs(r.x - ref.x);

        fdel << std::setprecision(16) << Delta << ","
             << EpsFix << ","
             << r.iters << ","
             << r.x << ","
             << err << "\n";

        std::cout << std::setprecision(6) << Delta << "\t"
                  << std::setprecision(0) << r.iters << "\t"
                  << std::setprecision(12) << r.x << "\t"
                  << std::setprecision(12) << err << "\n";
    }
    fdel.close();
    std::cout << "\nCSV сохранён: sensitivity_vs_delta.csv\n";

    return 0;
}
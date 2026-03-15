#include "methods.h"
#include <cmath>
#include <limits>

// ============================================================
// Функция Round - округление числа до заданной точности delta
// ============================================================
double Round(double value, double delta) {
    if (delta <= 0.0) {
        return value;
    }
    return std::round(value / delta) * delta;
}

// ============================================================
// Метод хорд
// Находит корень уравнения f(x) = 0 на отрезке [left, right]
// ============================================================
HordaResult HORDA(const std::function<double(double)>& f,
                  double left,
                  double right,
                  double eps,
                  int nmax,
                  double delta) {
    // Инициализация результата по умолчанию.
    HordaResult result{};
    result.x = std::numeric_limits<double>::quiet_NaN();
    result.a = left;
    result.b = right;
    result.iters = 0;
    result.ok = false;

    double a = left;
    double b = right;

    // Значения функции на концах отрезка с учетом возможного округления.
    double fa = Round(f(a), delta);
    double fb = Round(f(b), delta);

    // Проверка условия изоляции корня.
    if (!(fa * fb < 0.0)) {
        return result;
    }

    // Основной цикл метода хорд.
    for (int iter = 1; iter <= nmax; ++iter) {
        const double denominator = fb - fa;
        // Если знаменатель обнулился, следующую точку построить нельзя.
        if (denominator == 0.0) {
            result.x = 0.5 * (a + b);
            result.a = a;
            result.b = b;
            result.iters = iter;
            result.ok = false;
            return result;
        }

        // Точка пересечения хорды с осью Ox.
        const double c = a - fa * (b - a) / denominator;
        const double fc = Round(f(c), delta);

        // Остановка по малому значению функции или по длине отрезка.
        if (std::abs(fc) < eps || std::abs(b - a) < 2.0 * eps) {
            result.x = c;
            result.a = a;
            result.b = b;
            result.iters = iter;
            result.ok = true;
            return result;
        }

        // Оставляем ту часть отрезка, где сохраняется смена знака.
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // Если достигнут предел итераций, возвращаем последнее приближение.
    result.x = 0.5 * (a + b);
    result.a = a;
    result.b = b;
    result.iters = nmax;
    result.ok = false;
    return result;
}

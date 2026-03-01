#include "methods.h"
#include <cmath>
#include <limits>

double Round(double value, double Delta) {
    if (Delta <= 0.0) return value;               // без округления
    return std::round(value / Delta) * Delta;     // к ближайшему кратному Delta
}

BisectResult BISECT(const std::function<double(double)>& f,
                    double Left, double Right,
                    double Eps,
                    int nmax,
                    double Delta) {
    BisectResult res{};
    res.a = Left;
    res.b = Right;
    res.iters = 0;
    res.ok = false;
    res.x = std::numeric_limits<double>::quiet_NaN();

    double a = Left, b = Right;

    double fa = Round(f(a), Delta);
    double fb = Round(f(b), Delta);

    // Проверка условия теоремы Коши: f(a)*f(b) < 0
    if (!(fa * fb < 0.0)) {
        return res; // ok=false
    }

    for (int n = 1; n <= nmax; ++n) {
        double x = 0.5 * (a + b);
        double fx = Round(f(x), Delta);

        // Условие остановки по длине интервала: |b-a| < 2*Eps
        if (std::abs(b - a) < 2.0 * Eps) {
            res.x = x; res.a = a; res.b = b; res.iters = n; res.ok = true;
            return res;
        }

        // Если попали ровно в ноль (в т.ч. из-за округления)
        if (fx == 0.0) {
            res.x = x; res.a = a; res.b = b; res.iters = n; res.ok = true;
            return res;
        }

        // Выбираем половину, где знак меняется
        if (fa * fx < 0.0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }
    }

    // Не сошлось за nmax
    res.x = 0.5 * (a + b);
    res.a = a; res.b = b; res.iters = nmax; res.ok = false;
    return res;
}
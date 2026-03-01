#include "methods.h"
#include <cmath>
#include <limits>

double Round(double value, double Delta) {
    if (Delta <= 0.0) return value;
    return std::round(value / Delta) * Delta;
}

HordaResult HORDA(const std::function<double(double)>& f,
                  double Left, double Right,
                  double Eps,
                  int nmax,
                  double Delta) {
    HordaResult res{};
    res.a = Left;
    res.b = Right;
    res.iters = 0;
    res.ok = false;
    res.x = std::numeric_limits<double>::quiet_NaN();

    double a = Left, b = Right;
    double fa = Round(f(a), Delta);
    double fb = Round(f(b), Delta);

    if (!(fa * fb < 0.0)) {
        return res; // нет смены знака
    }

    for (int n = 1; n <= nmax; ++n) {
        // формула точки пересечения хорды с осью X:
        // c = a - f(a)*(b-a)/(f(b)-f(a))
        double denom = (fb - fa);
        if (denom == 0.0) { // защита
            res.x = 0.5 * (a + b);
            res.a = a; res.b = b; res.iters = n; res.ok = false;
            return res;
        }

        double c = a - fa * (b - a) / denom;
        double fc = Round(f(c), Delta);

        // условия остановки
        if (std::abs(fc) < Eps) {
            res.x = c; res.a = a; res.b = b; res.iters = n; res.ok = true;
            return res;
        }
        if (std::abs(b - a) < 2.0 * Eps) {
            res.x = c; res.a = a; res.b = b; res.iters = n; res.ok = true;
            return res;
        }

        // сохраняем отрезок с корнем (bracketing)
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    res.x = 0.5 * (a + b);
    res.a = a; res.b = b; res.iters = nmax; res.ok = false;
    return res;
}
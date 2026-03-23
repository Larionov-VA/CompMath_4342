#include "methods.h"
#include <cmath>
#include <limits>
#include <algorithm>

// Округление к ближайшему кратному Delta
double Round(double value, double Delta) {
    if (Delta <= 0.0) return value;
    return std::round(value / Delta) * Delta;
}

// Численная вторая производная (центральная разность)
// Нужна для оценки M2 и для выбора x0 (условие f(x0)*f''(x0)>0)
static double second_derivative_num(const std::function<double(double)>& f, double x) {
    double h = 1e-5 * std::max(1.0, std::abs(x));
    double f1 = f(x - h);
    double f2 = f(x);
    double f3 = f(x + h);
    return (f3 - 2.0 * f2 + f1) / (h * h);
}

NewtonResult NEWTON(const std::function<double(double)>& f,
                    const std::function<double(double)>& fp,
                    double Left, double Right,
                    double x0,
                    double Eps,
                    int nmax,
                    double Delta) {
    NewtonResult res{};
    res.x = std::numeric_limits<double>::quiet_NaN();
    res.iters = 0;
    res.ok = false;
    res.last_step = std::numeric_limits<double>::infinity();

    // Оценка m1 = min |f'| и M2 = max |f''| на [Left;Right] для критерия (3.1)
    double m1 = std::numeric_limits<double>::infinity();
    double M2 = 0.0;

    const int SAMPLES = 1000;
    for (int i = 0; i <= SAMPLES; ++i) {
        double x = Left + (Right - Left) * (double)i / (double)SAMPLES;
        double d1 = std::abs(fp(x));
        double d2 = std::abs(second_derivative_num(f, x));
        if (d1 < m1) m1 = d1;
        if (d2 > M2) M2 = d2;
    }

    // eps0 = 2*m1/M2 * Eps (формула (3.1))
    double eps0 = Eps;
    if (m1 > 0.0 && std::isfinite(m1) && M2 > 0.0 && std::isfinite(M2)) {
        eps0 = (2.0 * m1 / M2) * Eps;
    }

    double x_prev = x0;
    double fx_prev = Round(f(x_prev), Delta);
    double fpx_prev = Round(fp(x_prev), Delta);

    if (std::abs(fpx_prev) < 1e-15) {
        return res; // производная почти ноль => деление опасно
    }

    for (int n = 1; n <= nmax; ++n) {
        double x_next = x_prev - fx_prev / fpx_prev;
        double step = std::abs(x_next - x_prev);

        double fx_next  = Round(f(x_next), Delta);
        double fpx_next = Round(fp(x_next), Delta);

        res.last_step = step;

        // Основной критерий из методички: |x_n - x_{n-1}| < eps0
        if (step < eps0) {
            res.x = x_next;
            res.iters = n;
            res.ok = true;
            return res;
        }

        // Дополнительный контроль (не мешает, но полезен при больших Delta)
        if (std::abs(fx_next) < Eps) {
            res.x = x_next;
            res.iters = n;
            res.ok = true;
            return res;
        }

        if (std::abs(fpx_next) < 1e-15) {
            // дальше делить нельзя
            res.x = x_next;
            res.iters = n;
            res.ok = false;
            return res;
        }

        x_prev = x_next;
        fx_prev = fx_next;
        fpx_prev = fpx_next;
    }

    res.x = x_prev;
    res.iters = nmax;
    res.ok = false;
    return res;
}
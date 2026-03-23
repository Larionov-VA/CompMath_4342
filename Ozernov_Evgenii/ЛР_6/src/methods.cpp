#include "methods.h"
#include <cmath>
#include <limits>
#include <algorithm>

double Round(double value, double Delta) {
    if (Delta <= 0.0) return value;
    return std::round(value / Delta) * Delta;
}

double estimate_q(const std::function<double(double)>& dphi,
                  double Left, double Right,
                  int samples,
                  double Delta) {
    double q = 0.0;
    for (int i = 0; i <= samples; ++i) {
        double x = Left + (Right - Left) * static_cast<double>(i) / static_cast<double>(samples);
        double val = Round(dphi(x), Delta);
        q = std::max(q, std::abs(val));
    }
    return q;
}

IterationResult ITER(const std::function<double(double)>& phi,
                     const std::function<double(double)>& dphi,
                     double Left, double Right,
                     double x0,
                     double Eps,
                     int nmax,
                     double Delta) {
    IterationResult res{};
    res.x = std::numeric_limits<double>::quiet_NaN();
    res.iters = 0;
    res.ok = false;
    res.last_step = std::numeric_limits<double>::infinity();

    // Проверка условия сходимости q < 1
    double q = estimate_q(dphi, Left, Right, 2000, Delta);
    if (!(q < 1.0)) {
        return res;
    }

    double x_prev = Round(x0, Delta);

    for (int n = 1; n <= nmax; ++n) {
        double x_next = Round(phi(x_prev), Delta);
        double step = std::abs(x_next - x_prev);

        res.last_step = step;

        // Основной критерий из методички: |x_n - x_(n-1)| < Eps
        if (step < Eps) {
            res.x = x_next;
            res.iters = n;
            res.ok = true;
            return res;
        }

        // Небольшая защита: если ушли за интервал, считаем, что что-то пошло не так
        if (x_next < Left - 1e-12 || x_next > Right + 1e-12) {
            res.x = x_next;
            res.iters = n;
            res.ok = false;
            return res;
        }

        x_prev = x_next;
    }

    res.x = x_prev;
    res.iters = nmax;
    res.ok = false;
    return res;
}

#include "methods.h"

#include <cmath>

// ============================================================
// Округление значения до ближайшего числа, кратного delta
// ============================================================
double Round(double value, double delta) {
    if (delta <= 0.0) {
        return value;
    }
    return std::round(value / delta) * delta;
}

// ============================================================
// Оценка q = max |phi'(x)| на [left, right]
// ============================================================
double EstimateQ(const std::function<double(double)>& dphi,
                 double left,
                 double right,
                 int samples) {
    if (samples < 1) {
        samples = 1;
    }

    double q = 0.0;
    for (int i = 0; i <= samples; ++i) {
        const double t = static_cast<double>(i) / samples;
        const double x = left + (right - left) * t;
        const double value = std::fabs(dphi(x));
        if (value > q) {
            q = value;
        }
    }
    return q;
}

// ============================================================
// Метод простых итераций
// ============================================================
IterResult ITER(const std::function<double(double)>& phi,
                double x0,
                double eps,
                int nmax) {
    IterResult result{};
    result.x = x0;
    result.lastDiff = 0.0;
    result.iters = 0;
    result.ok = false;

    double x = x0;
    for (int n = 1; n <= nmax; ++n) {
        const double y = phi(x);
        const double diff = std::fabs(y - x);

        result.x = y;
        result.lastDiff = diff;
        result.iters = n;

        if (diff < eps) {
            result.ok = true;
            return result;
        }

        x = y;
    }

    return result;
}

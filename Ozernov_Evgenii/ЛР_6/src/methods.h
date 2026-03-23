#pragma once
#include <functional>

struct IterationResult {
    double x;         // найденное приближение корня
    int iters;        // число итераций
    bool ok;          // сошлось ли
    double last_step; // |x_n - x_(n-1)|
};

double Round(double value, double Delta);

// Оценка q = max |phi'(x)| на [Left; Right]
double estimate_q(const std::function<double(double)>& dphi,
                  double Left, double Right,
                  int samples = 1000,
                  double Delta = 0.0);

// Метод простых итераций:
// phi  - функция phi(x)
// dphi - производная phi'(x)
// Left, Right - интервал, где изолирован корень
// x0 - начальное приближение (лежит на [Left;Right])
// Eps - требуемая точность
// nmax - ограничение итераций
// Delta - округление значений phi и phi' через Round (0 => без округления)
IterationResult ITER(const std::function<double(double)>& phi,
                     const std::function<double(double)>& dphi,
                     double Left, double Right,
                     double x0,
                     double Eps,
                     int nmax = 1000000,
                     double Delta = 0.0);

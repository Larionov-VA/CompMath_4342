#pragma once
#include <functional>

struct NewtonResult {
    double x;        // найденное приближение корня
    int iters;       // число итераций
    bool ok;         // сошлось ли
    double last_step;// |x_n - x_{n-1}|
};

double Round(double value, double Delta);

// Метод Ньютона:
// f  - функция
// fp - производная f'(x)
// Left, Right - интервал, где изолирован корень (для оценки m1 и M2)
// x0 - начальное приближение (выбирается так, чтобы f(x0)*f''(x0)>0)
// Eps - требуемая точность
// nmax - ограничение итераций
// Delta - округление значений f и f' через Round (0 => без округления)
NewtonResult NEWTON(const std::function<double(double)>& f,
                    const std::function<double(double)>& fp,
                    double Left, double Right,
                    double x0,
                    double Eps,
                    int nmax = 1000000,
                    double Delta = 0.0);
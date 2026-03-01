#pragma once
#include <functional>

struct BisectResult {
    double x;        // найденное приближение корня
    double a, b;     // итоговый отрезок
    int iters;       // число итераций
    bool ok;         // сошлось ли (true/false)
};

// Округление значения до "точности Delta" (к ближайшему кратному Delta).
double Round(double value, double Delta);

// Метод бисекции.
// f      - функция
// Left,Right - начальный отрезок (должен давать разные знаки на концах)
// Eps    - требуемая точность
// nmax   - ограничение на число итераций
// Delta  - точность округления значений f(...) (0 => без округления)
BisectResult BISECT(const std::function<double(double)>& f,
                    double Left, double Right,
                    double Eps,
                    int nmax = 1000000,
                    double Delta = 0.0);
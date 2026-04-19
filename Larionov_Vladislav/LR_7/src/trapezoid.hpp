#include "./defines.hpp"

#define TRAPEZOID 4

ld trapezoid(std::function<ld(ld)> function, ld a, ld b, int n) {
    ld h = (b - a) / n;
    ld sum = 0.0L;
    ld fa = function(a);
    ld fb = function(b);
    if (fa == std::numeric_limits<ld>::max() || fb == std::numeric_limits<ld>::max()) {
        std::cerr << "Невозможно вычислить функцию на границах интервала.\n";
        exit(1);
    }
    sum = (fa + fb) / 2.0L;
    for (int i = 1; i < n; ++i) {
        ld x = a + i * h;
        ld fx = function(x);
        if (fx == std::numeric_limits<ld>::max()) {
            std::cerr << "Невозможно получить значение функции в точке " << x << '\n';
            exit(1);
        }
        sum += fx;
    }
    return h * sum;
}
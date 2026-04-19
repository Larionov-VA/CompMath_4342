#include "./defines.hpp"

#define SIMPSON 5


ld NewtonSimpson(std::function<ld(ld)> function, ld a, ld b, int n) {
    if (n % 2 != 0) {
        std::cerr << "Ошибка: для метода Симпсона число разбиений n должно быть чётным.\n";
        exit(1);
    }
    ld h = (b - a) / n;
    ld sum = 0.0L;
    ld fa = function(a);
    ld fb = function(b);
    if (fa == std::numeric_limits<ld>::max() || fb == std::numeric_limits<ld>::max()) {
        std::cerr << "Невозможно вычислить функцию на границах интервала.\n";
        exit(1);
    }
    sum = fa + fb;
    for (int i = 1; i < n; i += 2) {
        ld x = a + i * h;
        ld fx = function(x);
        if (fx == std::numeric_limits<ld>::max()) {
            std::cerr << "Невозможно получить значение функции в точке " << x << '\n';
            exit(1);
        }
        sum += 4.0L * fx;
    }
    for (int i = 2; i < n; i += 2) {
        ld x = a + i * h;
        ld fx = function(x);
        if (fx == std::numeric_limits<ld>::max()) {
            std::cerr << "Невозможно получить значение функции в точке " << x << '\n';
            exit(1);
        }
        sum += 2.0L * fx;
    }
    return (h / 3.0L) * sum;
}
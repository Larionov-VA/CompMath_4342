#include "./defines.hpp"

#define RECTANGLE 1

enum rectangleType {
    LEFT,
    RIGHT,
    CENTER
};

rectangleType rectMethodType;

ld rectangle(std::function<ld(ld)> function, ld a, ld b, int n) {
    ld result = 0.0L;
    ld h = (b - a) / (ld)n;
    int counter = 0;
    while (counter < n) {
        ld currentValue;
        ld point;
        switch (rectMethodType) {
        case rectangleType::RIGHT:
            point = a + h * ++counter;
            break;
        case rectangleType::CENTER:
            point = a + h * counter++ + h / 2;
            break;
        case rectangleType::LEFT:
            point = a + h * counter++;
            break;
        default:
            point = a + h * counter++ + h / 2;
        }
        currentValue = function(point);
        if (currentValue != std::numeric_limits<ld>::max()) {
            result += currentValue;
        }
        else {
            std::cerr << "Невозможно получить значение функции в точке "
            << std::to_string(point) << '\n';
            exit(1);
        }
    }
    return h*result;
}
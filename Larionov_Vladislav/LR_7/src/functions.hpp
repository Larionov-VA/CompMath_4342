#include "./defines.hpp"

#ifndef EPSILON
    #define EPSILON 1e-10
#endif

#define LINEAR 1
#define QUAD 2
#define POW 3
#define EXP 4
#define LOG 5
#define SIN 6
#define COS 7
#define TG 8
#define CTG 9


struct coefficients {
    ld a = 1, b, c;
};


coefficients funcCoeff;


ld linear(ld x) {
    return funcCoeff.a * x + funcCoeff.b;
}


ld quadratic(ld x) {
    return \
    funcCoeff.a * x * x +
    funcCoeff.b * x +
    funcCoeff.c;
}


ld power(ld x) {
    if (funcCoeff.b < 0 && std::abs(x) < EPSILON) {
        return std::numeric_limits<ld>::max();
    }
    else {
        return funcCoeff.a * std::pow(x, funcCoeff.b);
    }
}


ld exponential(ld x) {
    return funcCoeff.a * std::pow(funcCoeff.b, x);
}


ld lg(ld x) {
    if (funcCoeff.b + x <= 0) {
        return std::numeric_limits<ld>::max();
    }
    else {
        return funcCoeff.a + std::log2l(funcCoeff.b + x);
    }
}


ld sn(ld x) {
    return funcCoeff.a + std::sin(funcCoeff.b + x);
}


ld cs(ld x) {
    return funcCoeff.a + std::cos(funcCoeff.b + x);
}


ld tg(ld x) {
    ld arg = funcCoeff.b + x;
    if (std::abs(std::cos(arg)) < EPSILON) {
        return std::numeric_limits<ld>::max();
    }
    return funcCoeff.a + std::tan(arg);
}


ld cg(ld x) {
    ld arg = funcCoeff.b + x;
    if (std::abs(std::sin(arg)) < EPSILON) {
        return std::numeric_limits<ld>::max();
    }
    else {
        return funcCoeff.a + std::cos(arg) / std::sin(arg);
    }
}


ld mod2(ld x) {
    return std::abs(x - 1);
}
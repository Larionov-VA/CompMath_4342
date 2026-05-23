#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"
#define __NEWTON

double delta;

double F(double x) {
    double val = 5.0 * sin(x);
    if (val < 0) { puts("Ошибка: значение под корнем меньше нуля\n"); exit(1); }
    double res = sqrt(val);
    return Round(res, delta);
}

double derivativeF(double x) {
    double res = 2.0 * x - 5.0 * cos(x);
    return Round(res, delta);
}

int main(){
    double x0, eps, x;
    int k;
    float x0_in, eps1, delta1;

    printf("Лабораторная работа №6. Вариант 1\n");

    printf("Введите начальное приближение x0: "); // 2.1
    scanf("%f", &x0_in);
    x0 = (double)2.1;

    printf("Введите точность вычисления корня (Eps): "); // 0.0001
    scanf("%f", &eps1);
    eps = (double)eps1;

    printf("Введите точность задания функции (Delta): "); // 0.00000001
    scanf("%f", &delta1);
    delta = (double)delta1;

    x = ITER(x0, eps, k);

    printf("\nРезультаты метода итераций:\n");
    printf("Корень x = %7.7f\n", x);
    printf("Число итераций N = %d\n", k);

    return 0;

}

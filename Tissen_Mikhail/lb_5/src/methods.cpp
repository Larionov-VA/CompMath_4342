#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "methods.h"

static double gDelta = 0.0;
static double gM1 = 0.0;
static double gM2 = 0.0;

// Установка шага округления значений функции и производной.
void SetRoundDelta(double Delta)
{
    gDelta = Delta;
}

// Установка констант m1 и M2 для критерия остановки.
void SetNewtonConstants(double m1, double M2)
{
    gM1 = m1;
    gM2 = M2;
}

// ============================================================
// Функция Round - округление числа до заданной точности Delta
// ============================================================
double Round(double X, double Delta)
{
    if (Delta <= 0.0) return X;
    if (Delta < 1E-12) {
        printf("Неверное задание точности округления\n");
        exit(1);
    }

    if (X >= 0.0)
        return Delta * floor(X / Delta + 0.5);
    else
        return Delta * ceil(X / Delta - 0.5);
}

// ============================================================
// Метод Ньютона
// x_(n+1) = x_n - f(x_n) / f'(x_n)
// Критерий остановки согласно оценке погрешности:
// |x_n - x_(n-1)| < sqrt(2*m1*Eps/M2)
// ============================================================
double NEWTON(double X, double Eps, int &N)
{
    // Проверка корректности входных данных.
    if (Eps <= 0.0) {
        printf("Неверное задание точности\n");
        exit(1);
    }
    if (gM1 <= 0.0 || gM2 <= 0.0) {
        printf("Не заданы константы m1 и M2 для метода Ньютона\n");
        exit(1);
    }

    // Порог остановки по априорной оценке погрешности метода Ньютона.
    double E0 = sqrt(2.0 * gM1 * fabs(Eps) / gM2);
    double XPrev = X;
    double FX, DFX, XNext;

    N = 0;

    // Основной итерационный процесс метода Ньютона.
    while (true) {
        FX = Round(F(XPrev), gDelta);
        DFX = Round(DF(XPrev), gDelta);

        // Производная в текущей точке не должна обращаться в нуль.
        if (DFX == 0.0) {
            printf("Производная равна нулю. Метод Ньютона неприменим.\n");
            exit(1);
        }

        XNext = XPrev - FX / DFX;
        N++;

        // Критерий остановки по разности соседних приближений.
        if (fabs(XNext - XPrev) < E0)
            return XNext;

        // Новое приближение должно оставаться в области определения.
        if (XNext <= 0.0) {
            printf("Следующее приближение вышло из области определения\n");
            exit(1);
        }

        // Защита от зацикливания при нарушении сходимости.
        if (N > 1000) {
            printf("Превышено допустимое число итераций\n");
            exit(1);
        }

        XPrev = XNext;
    }
}

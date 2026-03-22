#include "methods.hpp"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef ROUND
    #define ROUND 0
#endif
#ifndef DELTA
    #define DELTA 0.001
#endif

#if !ROUND
double F(double x) {
    return std::pow(x, 3) - 3 * x - 2 * std::pow(M_El, -x);
}
double derivativeF(double x) {
    return 2 * std::pow(M_El, -x) + 3 * std::pow(x, 2) - 3;
}
double phi(double x) {
    return pow(3*x + 2*exp(-x), 1.0/3.0);
}
double dphi(double x) {
    double temp = 3*x + 2*exp(-x);
    double denom = 3 * pow(temp, 2.0/3.0);
    return (3 - 2*exp(-x)) / denom;
}
#else
double F(double x) {
    return Round(std::pow(x, 3) - 3 * x - 2 * std::pow(M_El, -x), DELTA);
}
double derivativeF(double x) {
    return Round(2 * std::pow(M_El, -x) + 3 * std::pow(x, 2) - 3, DELTA);
}
double phi(double x) {
    return Round(pow(3*x + 2*exp(-x), 1.0/3.0), DELTA);
}
double dphi(double x) {
    double temp = 3*x + 2*exp(-x);
    double denom = 3 * pow(temp, 2.0/3.0);
    return Round((3 - 2*exp(-x)) / denom, DELTA);
}
#endif

double BISECT(double Left, double Right, double Eps, int &N)
{
    double E = fabs(Eps) * 2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left + Right) / 2.0;
    double Y;
    if (FLeft * FRight > 0.0)
    {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0)
    {
        puts("Неверное задание точности\n");
        exit(1);
    }
    N = 0;
    if (FLeft == 0.0)
        return Left;
    if (FRight == 0.0)
        return Right;
    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left);
        /* вычисление середины отрезка
         */
        Y = F(X);
        if (Y == 0.0) {
            if (derivativeF(X)) {
                return (X);
            }
            else {
                puts("Производная не равна нулю\n");
                exit(1);
            }
        }
        if (Y * FLeft < 0.0)
            Right = X;
        else
        {
            Left = X;
            FLeft = Y;
        }
        N++;
    };
    if (!derivativeF(X)) {
        puts("Производная не равна нулю\n");
        exit(1);
    }
    return (X);
}

double Round(double X, double Delta)
{
    if (Delta <= 1E-9)
    {
        puts("Неверное задание точности округления\n");
        exit(1);
    }
    if (X > 0.0)
        return (Delta * (long((X / Delta) + 0.5)));
    else
        return (Delta * (long((X / Delta) - 0.5)));
}

double ITER(double X0, double Eps, int &N)
{
    if (Eps <= 0.0)
    {
        puts("Неверное задание точности\n");
        exit(1);
    }
    double X1 = phi(X0);
    double X2 = phi(X1);
    N = 2;
    while ((X1 - X2) * (X1 - X2) > fabs((2 * X1 - X0 - X2) * Eps))
    {
        X0 = X1;
        X1 = X2;
        X2 = phi(X1);
        N++;
    }
    return (X2);
}

#ifdef __NEWTON
double NEWTON(double X, double Eps, int &N)
{
    double Y, Y1, DX;
    N = 0;
    do
    {
        Y = F(X);
        if (Y == 0.0)
            return (X);
        Y1 = derivativeF(X);
        if (Y1 == 0.0)
        {
            puts("Производная обратилась в ноль\n");
            exit(1);
        }
        DX = Y / Y1;
        X = X - DX;
        N++;
    } while (fabs(DX) > Eps);
    return (X);

}
#endif

double HORDA(double Left, double Right, double Eps, int &N)
{
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;
    if (FLeft * FRight > 0.0)
    {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0)
    {
        puts("Неверное задание точности\n");
        exit(1);
    }
    N = 0;
    if (FLeft == 0.0)
        return Left;
    if (FRight == 0.0)
        return Right;
    do
    {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        if (Y == 0.0)
            return (X);
        if (Y * FLeft < 0.0)
        {
            Right = X;
            FRight = Y;
        }
        else
        {
            Left = X;
            FLeft = Y;
        }
        N++;
    } while (fabs(Y) >= Eps);
    return (X);
}
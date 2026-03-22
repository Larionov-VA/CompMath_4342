#include "methods.hpp"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef ROUND_FLAG
    #define ROUND_FLAG 0
#endif
#ifndef DELTA
    #define DELTA 0.001
#endif

#if !ROUND_FLAG
extern double F(double x) {
    // Для метода простых итераций (функция ITER) F(x) должна возвращать эквивалентную функцию Phi(x), а не саму f(x)=0.
    // Преобразование: x^2 = arctg(1/x) => x = sqrt(arctg(1/x))
    return sqrt(atan(1.0 / x));
}
#ifdef __NEWTON
extern double F1(double x) {
    return (-1.0 / (x * x + 1.0)) - (2.0 * x);
}
#endif
#else
extern double F(double x) {
    return Round(sqrt(atan(1.0 / x)), DELTA);
}
#ifdef __NEWTON
extern double F1(double x) {
    return Round((-1.0 / (x * x + 1.0)) - (2.0 * x), DELTA);
}
#endif
#endif

double HORDA(double Left, double Right, double Eps, int &N)
{
    double FLeft = F(Left);
    double FRight = F(Right);
    double X, Y;
    
    if (FLeft * FRight > 0.0) {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0) {
        puts("Неверное задание точности\n");
        exit(1);
    }
    
    N = 0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    
    do
    {
        X = Left - (Right - Left) * FLeft / (FRight - FLeft);
        Y = F(X);
        if (Y == 0.0) return (X);
        
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
    }
    while ( fabs(Y) >= Eps );
    
    return(X);
}

double BISECT(double Left, double Right, double Eps, int &N)
{
    double E = fabs(Eps) * 2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left + Right) / 2.0;
    double Y;
    if (FLeft * FRight > 0.0) {puts("Неверное задание интервала\n");exit(1);}
    if (Eps <= 0.0) {puts("Неверное задание точности\n");exit(1);}
    N=0;
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left);
        Y = F(X);
        if (Y == 0.0) return (X);
        if (Y * FLeft < 0.0) Right = X;
        else { Left = X; FLeft = Y; }
        N++;
    };
    return(X);
}

double Round(double X, double Delta)
{
    if (Delta <= 1E-9) {
        puts("Неверное задание точности округления\n");
        exit(1);
    }
    if (X > 0.0) return (Delta * (long((X / Delta) + 0.5)));
    else         return (Delta * (long((X / Delta) - 0.5)));
}

double ITER(double X0, double Eps, int &N)
{
    if (Eps <= 0.0) {puts("Неверное задание точности\n"); exit (1);}
    // Здесь F() работает как Phi(), возвращая новое X
    double X1 = F(X0);
    double X2 = F(X1);
    N = 2;
    // Особенный критерий остановки, использующий 3 точки (процесс Эйткена)
    while( (X1 - X2)*(X1 - X2) > fabs((2*X1 - X0 - X2)*Eps) )
    {
        X0 = X1;
        X1 = X2;
        X2 = F(X1);
        N++;
    }
    return(X2);
}

#ifdef __NEWTON
double NEWTON (double X, double Eps, int &N)
{
    extern double F1 (double);
    double Y, Y1, DX;
    N = 0;
    do
    {
        Y = F(X);
        if (Y == 0.0) return (X);
        
        Y1 = F1(X);
        if (Y1 == 0.0) {
            puts("Производная обратилась в ноль\n"); 
            exit(1);
        }
        
        DX = Y / Y1; 
        X = X - DX; 
        N++;
    }
    while (fabs(DX) > Eps);
    
    return (X);
}
#endif
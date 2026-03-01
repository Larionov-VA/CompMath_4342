/*********************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern double F(double);
/*
МЕТОД БИСЕКЦИИ (Половинного деления)
Left, Right - границы интервала
Eps - точность
N - переменная для возврата количества итераций
*/
double BISECT(double Left, double Right, double Eps, int &N)
{
    double E = fabs(Eps) * 2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left + Right) / 2.0;
    double Y;

    if (FLeft * FRight > 0.0) {
        puts("Неверное задание интервала\n");
        exit(1);
    }
    if (Eps <= 0.0) {
        puts("Неверное задание точности\n");
        exit(1);
    }

    N = 0;
    if (FLeft == 0.0)  return Left;
    if (FRight == 0.0) return Right;

    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left); // Вычисление середины отрезка
        Y = F(X);

        if (Y == 0.0) return (X);

        if (Y * FLeft < 0.0) {
            Right = X;
        } else {
            Left = X;
            FLeft = Y;
        }
        N++;
    }
    return (X);
}

/*
ФУНКЦИЯ ОКРУГЛЕНИЯ
X - число
Delta - шаг округления (например, 0.01)
*/
double Round(double X, double Delta)
{
    if (Delta <= 1E-9) {
        puts("Неверное задание точности округления\n");
        exit(1);
    }
    if (X > 0.0) 
        return (Delta * (long((X / Delta) + 0.5)));
    else         
        return (Delta * (long((X / Delta) - 0.5)));
}

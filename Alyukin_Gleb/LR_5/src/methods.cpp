#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern double F(double);
extern double derivativeF(double);

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
        DX = Y / Y1; // DX - Разность двух последних шагов
        X = X - DX; // DX = X_старое - X_новое
        N++;
    } while (fabs(DX) > Eps);
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

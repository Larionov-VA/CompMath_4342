#include <stdio.h>
#include <math.h>
#include <stdlib.h>

extern double F(double);
extern double derivativeF(double);

double ITER(double X0,double Eps,int &N)
{
    if (Eps<=0.0) {puts("Неверное задание точности\n");exit (1);}
    double X1=F(X0);
    double X2=F(X1);
    N = 2;
    // разность соседних приближений |xn-1 - xn| < Eps
    while( (X1 - X2)*(X1 - X2) > fabs((2*X1-X0-X2)*Eps) && (N < 500))
    {
        X0 = X1;
        X1 = X2;
        X2 = F(X1);
        N++;

        if (X1 == X2) break;
    }
    return(X2);
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


double CurrentDelta = 0.0;

double F(double x);
double Round(double X, double Delta);
double BISECT(double Left, double Right, double Eps, int &N);

double F(double x)
{
    double raw = exp(-x) - x * x * x;  
    if (CurrentDelta > 0.0)
        return Round(raw, CurrentDelta);  
    else
        return raw;
}

double Round(double X, double Delta)
{
    if (Delta <= 1E-9)
    {
        puts("Неверное задание точности округления\n");
        exit(1);
    }
    if (X > 0.0)
        return Delta * (long)((X / Delta) + 0.5);
    else
        return Delta * (long)((X / Delta) - 0.5);
}

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
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;

    while ((Right - Left) >= E)
    {
        X = 0.5 * (Right + Left);        
        Y = F(X);

        if (Y == 0.0) return X;           

        if (Y * FLeft < 0.0)              
            Right = X;
        else                               
        {
            Left = X;
            FLeft = Y;                    
        }
        N++;
    }
    return X;                            
}

int main()
{
    double Left = 0.5;
    double Right = 1.0;

    double EpsValues[] = {0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001};
    int numEps = sizeof(EpsValues) / sizeof(EpsValues[0]);

    printf("%d\n", numEps);

    CurrentDelta = 0.0;   

    for (int i = 0; i < numEps; i++)
    {
        double eps = EpsValues[i];
        int iter;
        double x = BISECT(Left, Right, eps, iter);
        printf("%.15lf %d %.15lf\n", eps, iter, x);
    }

    double EpsFixed = 1e-6;
    double DeltaValues[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    int numDelta = sizeof(DeltaValues) / sizeof(DeltaValues[0]);

    printf("%d\n", numDelta);

    for (int i = 0; i < numDelta; i++)
    {
        CurrentDelta = DeltaValues[i];
        int iter;
        double root = BISECT(Left, Right, EpsFixed, iter);
        printf("%.15lf %d %.15lf\n", CurrentDelta, iter, root);
    }

    return 0;
}
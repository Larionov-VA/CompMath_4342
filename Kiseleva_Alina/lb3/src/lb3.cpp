#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>


double F(double x){
    return tan(x) - 1/x;
}

double BISECT(double Left, double Right, double Eps, int &N){
    double E = std::abs(Eps)*2.0;
    double FLeft = F(Left);
    double FRight = F(Right);
    double X = (Left+Right)/2.0;
    double Y;
    
    if (FLeft*FRight>0.0) {
        std::cout << "Неверное задание интервала\n";
        exit(1);
    }
    
    if (Eps<=0.0) {
        std::cout << "Неверное задание точности\n";
        exit(1);
    }

    N=0;
    
    if (FLeft==0.0) return Left;
    if (FRight==0.0) return Right;
    
    while ((Right-Left)>=E){
        X = 0.5*(Right + Left);
        Y = F(X);
        
        if (Y == 0.0) return (X);
        
        if (Y*FLeft < 0.0)
            Right=X;
        else { 
            Left=X;
            FLeft=Y; 
        }
        N++;
    }
    return(X);
}


double Round (double X, double Delta){
    if (Delta<=1E-9) {
        std::cout << "Неверное задание точности округления\n";
        exit(1);
    }
    
    if (X>0.0) 
        return (Delta*(long((X/Delta)+0.5)));
    else
        return (Delta*(long((X/Delta)-0.5)));
}

int main(){
    int N;
    double left = 0.5;
    double right = 1.0;

    std::ofstream file("results.txt");
    file << "Количество итераций\tТочность вычисления\tТочность округления\tКорень\n";

    for (double e = 0.1; e >= 0.000001; e *= 0.1){
        int cntSignesE = static_cast<int>(-log10(e));

        double root = BISECT(left, right, e, N);
        file << std::fixed << std::setprecision(cntSignesE);
        file << N << '\t' << e << '\t' << 0 << '\t' << root << '\n';
        
        for (double delta = 0.1; delta > 1E-9; delta *= 0.1){
            double rootWithRound = BISECT(Round(left, delta), Round(right, delta), e, N);

            int cntSignesDelta = static_cast<int>(-log10(delta));
            file << std::fixed << std::setprecision(cntSignesE);
            file << N << '\t' << e << '\t';
            file << std::fixed << std::setprecision(cntSignesDelta);
            file << delta << '\t';
            file << std::fixed << std::setprecision(cntSignesE+cntSignesDelta);
            file << rootWithRound << '\n';
        }
    }

    file.close();

    return 0;
}
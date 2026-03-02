#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>


double F(double x){
    return tan(x) - 1/x;
}

double HORDA(double Left, double Right, double Eps, std::vector<double>& N){
    double FLeft = F(Left);
    double FRight = F(Right);
    double X,Y;

    if (FLeft*FRight>0.0) {
        std::cout << "Неверное задание интервала\n";
        exit(1);
    }
    
    if (Eps<=0.0) {
        std::cout << "Неверное задание точности\n";
        exit(1);
    }
    
    if (FLeft==0.0) return Left;
    if (FRight==0.0) return Right;
    
    double x0 = 0.860334;
    std::vector<double> roots;

    do{
        X = Left - (Right - Left)*FLeft / (FRight - FLeft);
        Y = F(X);
        
        if (Y == 0.0) return (X);
        
        if (Y*FLeft < 0.0){ 
            Right=X; 
            FRight=Y; 
        } else { 
            Left=X; 
            FLeft=Y; 
        }
        
        double c = 0.0;

        if (roots.size()){
            c = fabs(X - x0) / fabs(roots.back() - x0);
        }

        roots.push_back(X);
        N.push_back(c);

    } while ( fabs(Y) >= Eps );
    
    return(X);
}

int main(){
    double left = 0.5;
    double right = 1.0;
    std::vector<double> N;

    std::ofstream f("convergenceSpeed.txt");
    f << "№ итерации\tСкорость сходимости\n";
    f << std::fixed << std::setprecision(6);

    double root = HORDA(left, right, 1e-5, N);

    for (int i = 0; i < N.size(); i++){
        f << i+1 << '\t' << N[i] << '\n';
    }

    f.close();

    return 0;
}
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

double f(double x){
    return tan(x) - 1/x;
}

double f1(double x){    //Производная
    return 1/(cos(x)*cos(x)) + 1/(x*x);
}

double F(double x){     //Фи для простых итераций
    return x - 0.855*f(x);
}

double bisect(double Left, double Right, double Eps, int &N){
    double E = std::abs(Eps)*2.0;
    double FLeft = f(Left);
    double FRight = f(Right);
    double X = (Left + Right) / 2.0;
    double Y;
    
    if (FLeft * FRight > 0.0) {
        std::cout << "Неверное задание интервала\n";
        exit(1);
    }
    
    if (Eps <= 0.0) {
        std::cout << "Неверное задание точности\n";
        exit(1);
    }

    N = 0;
    
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    
    while ((Right - Left) >= E){
        X = 0.5*(Right + Left);
        Y = f(X);
        
        if (Y == 0.0) return (X);
        
        if (Y*FLeft < 0.0)
            Right = X;
        else { 
            Left = X;
            FLeft = Y; 
        }
        
        N++;
    }
    
    return(X);
}

double horda(double Left, double Right, double Eps, int &N){
    double FLeft = f(Left);
    double FRight = f(Right);
    double X,Y;

    if (FLeft * FRight > 0.0) {
        std::cout << "Неверное задание интервала\n";
        exit(1);
    }
    
    if (Eps <= 0.0) {
        std::cout << "Неверное задание точности\n";
        exit(1);
    }
    
    N = 0;
    
    if (FLeft == 0.0) return Left;
    if (FRight == 0.0) return Right;
    
    do{
        X = Left - (Right - Left)*FLeft / (FRight - FLeft);
        Y = f(X);
        
        if (Y == 0.0) return (X);
        
        if (Y*FLeft < 0.0){ 
            Right = X; 
            FRight = Y; 
        } else { 
            Left = X; 
            FLeft = Y; 
        }
        N++;

    } while ( fabs(Y) >= Eps );
    
    return(X);
}

double newton(double x, double eps, int &cntIter){
    double y, y1, dx;
    cntIter = 0;

    do {
        y = f(x);
        if (y == 0.0) 
            return (x);

        y1 = f1(x);
        
        if (y1 == 0.0){
            std::cout << "Производная обратилась в ноль\n";
            exit(1);
        }
        
        dx = y / y1; 
        x = x - dx; 
        cntIter++;

    } while (fabs(dx) > eps);

    return x;
}

double iter(double x0, double eps, int& n){
    if (eps <= 0.0){
        std::cout << "Неверное задание точности\n";
        exit(1);
    }

    double x1 = F(x0);
    double x2 = F(x1);

    n = 2;

    while((x1 - x2)*(x1 - x2) > fabs((2*x1 - x0 - x2) * eps)){
        x0 = x1;
        x1 = x2;
        x2 = F(x1);
        n++;
    }

    return(x2);
}

double Round(double x, double delta){
    if (delta <= 1E-9) {
        std::cout << "Неверное задание точности округления\n";
        exit(1);
    }
    
    if (x > 0.0) 
        return (delta * (long((x/delta) + 0.5)));
    else
        return (delta * (long((x/delta) - 0.5)));
}

int main(){
    int N;
    double left = 3.3;  //Для бисекции и хорд
    double right = 3.5;
    double x0 = 3.5;    //Для Ньютона и простых итераций

    std::vector<std::string> methods = {
        "========== МЕТОД БИСЕКЦИИ ==========\n\n",
        "========== МЕТОД ХОРД ==========\n\n",
        "========== МЕТОД НЬЮТОНА ==========\n\n",
        "========== МЕТОД ПРОСТЫХ ИТЕРАЦИЙ ==========\n\n"
    };

    int switcher;
    std::cout << "Введите номер метода (0-3)\n";
    std::cin >> switcher;

    std::ofstream file("results.txt");
    file << methods[switcher];
    file << "Количество итераций\tТочность вычисления\tТочность округления\tКорень\n";

    for (double e = 0.1; e >= 1E-9; e *= 0.1){
        int cntSignesE = static_cast<int>(-log10(e));
        double root;

        file << std::fixed << std::setprecision(cntSignesE);

        if (switcher == 0){
            root = bisect(left, right, e, N);
            file << N << '\t' << e << '\t' << 0 << '\t' << root << '\n';
        
        } else if (switcher == 1){
            root = horda(left, right, e, N);
            file << N << '\t' << e << '\t' << 0 << '\t' << root << '\n';
        
        } else if (switcher == 2){
            root = newton(x0, e, N);
            file << N << '\t' << e << '\t' << 0 << '\t' << root << '\n';
        
        } else if (switcher == 3){
            root = iter(x0, e, N);
            file << N << '\t' << e << '\t' << 0 << '\t' << root << '\n';
        }
        
        for (double delta = 0.1; delta > 1E-9; delta *= 0.1){
            double rootWithRound = Round(root, delta); 

            int cntSignesDelta = static_cast<int>(-log10(delta));
            //file << std::fixed << std::setprecision(cntSignesE);
            //file << N << '\t' << e << '\t';
            //file << std::fixed << std::setprecision(cntSignesDelta);
            //file << delta << '\t';
            file << std::fixed << std::setprecision(cntSignesE+cntSignesDelta);
            file << rootWithRound << '\n';
        }
    }

    file.close();

    return 0;
}
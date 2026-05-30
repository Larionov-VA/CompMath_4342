#include <algorithm>
#include <iostream>
#include <string_view>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>

#define EPS 0.001


double calculateDistance(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return std::sqrt(sum);
}


void YACOBY(int n, std::vector<std::vector<double>>& A, 
        std::vector<double>& b, std::vector<double>& x, int &k){
    
    std::vector<double> x_new(n);

    while(1){
        
        for (int i = 0; i < n; i++){
            double S = 0.0;
            for (int j = 0; j < n; j++){
                if (i != j){
                    S += A[i][j] * x[j]; 
                }
            }

            x_new[i] = (b[i] - S) / A[i][i];
        }

        if (calculateDistance(x, x_new) < EPS){
            x = x_new;
            break;
        }

        x = x_new;
        k++;
    }

    return;
}


std::vector<std::vector<double>> generateDiagonallyDominantMatrix(int n, double min_val = -5.0, double max_val = 5.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                A[i][j] = dis(gen);
                sum += std::abs(A[i][j]);
            }
        }
        A[i][i] = sum + 1.0; 
    }
    return A;
}


std::vector<double> generateBFromReference(int n, const std::vector<std::vector<double>>& A) {
    std::vector<double> b(n, 0.0);
    std::vector<double> x_reference(n, 0.0);

    for (int j = 1; j < n; j++){
        x_reference[j] = x_reference[j-1] + 1.0;
    }

    std::cout << "Reference X: ";
    std::for_each(x_reference.begin(), x_reference.end(), [](double x){
                std::cout << x << ' ';
            });

    std::cout << std::endl;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += A[i][j] * x_reference[j];
        }
    }
    return b;
}

bool is_correct(std::vector<std::vector<double>>& A){
    for (int i = 0; i < A.size(); i++){
        if (A[i][i] == 0){
            return false;
        }
    }

    return true;
}


void tests(){
    std::ofstream out("jakoby.txt"); 
    if (!out.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
    }

    std::vector<int> size = {3, 10, 50, 100};

    for (int n : size){
        
        std::cout << "Size: " << n << std::endl;

        std::vector<double> x(n);
        std::vector<std::vector<double>> A = generateDiagonallyDominantMatrix(n);
    
        while (!is_correct(A)){
            A = generateDiagonallyDominantMatrix(n);
    
        }


        out << "Matrix A: " << std::endl;

        std::for_each(A.begin(), A.end(), [&](std::vector<double> x){
                    std::for_each(x.begin(), x.end(), [&](double a){
                                out << a << ' ';
                            });
                    out << std::endl;
                });

        out << std::endl;

        std::vector<double> b = generateBFromReference(n, A); 

        int k = 0;

        out << "Vector B: " ;
        std::for_each(b.begin(), b.end(), [&](double x){
                out << x << ' ';
            });

        out << std::endl << std::endl;

        YACOBY(n, A, b, x, k);


        // Print result
        std::cout << "Result: ";
        std::for_each(x.begin(), x.end(), [](double x){
                    std::cout << x << ' '; 
                });
        
        std::cout << std::endl;
        
        std::cout << "Error: ";
        
        
        std::vector<double> x_reference(n, 0.0);

        for (int j = 1; j < n; j++){
            x_reference[j] = x_reference[j-1] + 1.0;
        }

        std::cout << calculateDistance(x, x_reference) << std::endl;

        std::cout << std::endl << "Count iterration: ";
        std::cout << k << std::endl;

    }
}


int main(int argc, char* argv[]){

    for (int i = 0; i < argc; ++i) {
        if (std::string_view(argv[i]) == "-t"){
            tests();
            return 0;
        }

    }

    std::ofstream out("jakoby.txt"); 
    if (!out.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    int n;
    std::cin >> n;

    std::vector<double> x(n);
    std::vector<std::vector<double>> A = generateDiagonallyDominantMatrix(n);
    
    while (!is_correct(A)){
        A = generateDiagonallyDominantMatrix(n);
    
    }
    
    out << "Matrix A: " << std::endl;

    std::for_each(A.begin(), A.end(), [&](std::vector<double> x){
                std::for_each(x.begin(), x.end(), [&](double a){
                            out << a << ' ';
                        });
                out << std::endl;
            });

    out << std::endl;

    std::vector<double> b = generateBFromReference(n, A); 

    out << "Vector B: " ;
    std::for_each(b.begin(), b.end(), [&](double x){
                out << x << ' ';
            });

    out << std::endl << std::endl;
    
    int k = 0;


    YACOBY(n, A, b, x, k);


    // Print result
    std::cout << "Result: ";
    std::for_each(x.begin(), x.end(), [](double x){
                std::cout << x << ' '; 
            });

    std::cout << std::endl;

    
    std::cout << "Error: ";
        
        
    std::vector<double> x_reference(n, 0.0);

    for (int j = 1; j < n; j++){
        x_reference[j] = x_reference[j-1] + 1.0;
    }

    std::cout << calculateDistance(x, x_reference) << std::endl;

    std::cout << std::endl << "Count iterration: ";
    std::cout << k << std::endl;

    return 0;
}

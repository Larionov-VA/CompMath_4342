#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <lapacke.h>

int N;
const int n = 10000;
double matrix[n][n+1];
double eps = 1e-12;


void generateMatrix(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 10.0); //Равномерное распределение

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < (N+1); ++j){
            matrix[i][j] = dist(gen);
        }
    }
}

bool forwardPhase(){
    int currentRow = 0;     //Номер строки, где должен лежать следующий опорный элемент

    for (int col = 0; col < N; ++col){
        int pivotRow = currentRow;
        double maxVal = fabs(matrix[currentRow][col]);

        for (int row = currentRow + 1; row < N; ++ row){    //Поиск максимального элемента в столбце
            if (fabs(matrix[row][col]) > maxVal){
                maxVal = fabs(matrix[row][col]);
                pivotRow = row;
            }
        }

        if (maxVal < eps) return false;

        if (pivotRow != currentRow){    //Перестановка строк, если максимальный элемент не был в опорной строке
            for (int i = 0; i < N + 1; ++i){
                double temp = matrix[currentRow][i];
                matrix[currentRow][i] = matrix[pivotRow][i];
                matrix[pivotRow][i] = temp;
            }
        }

        for (int row = currentRow + 1; row < N; ++row){
            double mult = matrix[row][col]/matrix[currentRow][col];
            for (int j = col; j < N + 1; ++j){
                matrix[row][j] -= mult * matrix[currentRow][j];
            }
        }

        currentRow++;
    }

    return true;
}

std::pair<double, double> calcNormas(std::vector<double>& ansGauss, std::vector<double>& ansLapacke){
    double oneNorma = 0.0;
    double twoNorma = 0.0;
    double oneNormaLapacke = 0.0;
    double twoNormaLapacke = 0.0;

    for (int i = 0; i < N; ++i){
        oneNorma += fabs(ansLapacke[i] - ansGauss[i]);
        twoNorma += (ansLapacke[i] - ansGauss[i]) * (ansLapacke[i] - ansGauss[i]);
        oneNormaLapacke += fabs(ansLapacke[i]);
        twoNormaLapacke += ansLapacke[i] * ansLapacke[i];
    }

    return {oneNorma / oneNormaLapacke, sqrt(twoNorma) / sqrt(twoNormaLapacke)};
}

std::vector<double> backwardPhase(){
    std::vector<double> ans(N);

    for (int i = N-1; i >= 0; --i){
        double sum = 0.0;
        for (int j = i+1; j < N; ++j){
            sum += matrix[i][j] * ans[j];   //Вычисление суммы произведений коэффициентов и уже найденных X
        }
        ans[i] = (matrix[i][N] - sum) / matrix[i][i];
    }

    return ans;
}

bool gauss(std::vector<double>& ans){
    if (!forwardPhase()){
        std::cout << "Система не имеет решений или их бесконечное количество.\n";
        return false;
    }

    ans = backwardPhase();
    return true;
}


void lapack(std::vector<double>& A, std::vector<double>& b){
    int ipiv[N];
    int info = LAPACKE_dgesv(
            LAPACK_ROW_MAJOR,  // порядок хранения (по строкам)
            N,                 // размер матрицы
            1,                 // одна правая часть
            A.data(),          // матрица A
            N,                 // ведущая размерность
            ipiv,       // pivot-индексы
            b.data(),          // правая часть (на выходе здесь будет решение)
            1                  // ведущая размерность b 
        );
}

void compareGaussLapacke(std::vector<double>& matrixForLapacke, std::vector<double>& rightPart){
    std::vector<double> ansGauss(N);
    if (gauss(ansGauss)){
        lapack(matrixForLapacke, rightPart);
        
        std::pair<double, double> normas = calcNormas(ansGauss, rightPart);

        std::cout << "Относительная погрешность метода Гаусса по 1-норме = " << normas.first << '\n';
        std::cout << "Относительная погрешность метода Гаусса по 2-норме = " << normas.second << '\n';
        std::cout << "За точное взято решение с помощью Lapacke\n\n";
    }
}

int main(){
    std::cout << "КОМПИЛЯЦИЯ: g++ gauss.cpp -o gauss -llapacke -llapack -lblas\n";
    std::cout << "Введите размер матрицы <= 10 000\n";
    std::cin >> N;

    generateMatrix();

    std::vector<double> matrixForLapacke(N*N);
    std::vector<double> rightPart(N);

    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N+1; ++j){
            if (j == N) rightPart[i] = matrix[i][j];
            else matrixForLapacke[i*N+j] = matrix[i][j];
        }
    }

    std::cout << "Случайная матрица размера " << N << "\n\n";
   
    compareGaussLapacke(matrixForLapacke, rightPart);

    /*std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 10.0);

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N+1; ++j){
            if (j != N){
                matrix[i][j] = 1.0 / (i + j + 1);   //Матрица Гильберта
                matrixForLapacke[i*N+j] = matrix[i][j];
            } else {
                matrix[i][j] = dist(gen);
                rightPart[i] = matrix[i][j];
            }
        }
    }

    std::cout << "Матрица Гильберта размера " << N << "\n\n";
    
    compareGaussLapacke(matrixForLapacke, rightPart);*/

    return 0;
}
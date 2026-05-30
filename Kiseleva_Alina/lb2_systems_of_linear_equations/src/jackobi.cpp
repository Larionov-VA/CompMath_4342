#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <lapacke.h>

int N;
const int n = 15000;
double matrix[n][n+1];
double eps = 1e-2;


void generateMatrix(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 10.0); //Равномерное распределение
    std::uniform_real_distribution<double> distDiag(10.0*N, 11.0*N);

    for (int i = 0; i < N; ++i){
        for (int j = 0; j < (N+1); ++j){
            if (i == j) matrix[i][j] = distDiag(gen);
            else matrix[i][j] = dist(gen);
        }
    }
}

double normaForJakobi(std::vector<double>& ansPrev, std::vector<double>& ansCurr){
    double norma = 0.0;
    for (int i = 0; i < N; ++i){
        norma += (ansCurr[i] - ansPrev[i]) * (ansCurr[i] - ansPrev[i]);
    }
    return sqrt(norma);
}

bool checkDiagonalDominance(){
    for (int i = 0; i < N; ++i){
        double sum = 0.0;

        for (int j = 0; j < N; ++j){
            if (j != i)
                sum += fabs(matrix[i][j]);
        }

        if (fabs(matrix[i][i]) <= sum)
            return false;
    }

    return true;
}

bool jakobi(std::vector<double>& ans, int& cntIter){
    if (!checkDiagonalDominance()){
        std::cout << "Метод Якоби не применим для решения данной системы.\n";
        return false;

    } else {
        std::vector<double> xOld(N);   

        double residual = 10.0;
        cntIter = 0;

        while (residual > eps){
            std::vector<double> xCurr(N);

            for (int i = 0; i < N; ++i){
                double b = matrix[i][N];
                double sum = 0.0;

                for (int j = 0; j < N; ++j){
                    if (i != j) 
                        sum += matrix[i][j] * xOld[j];
                }

                xCurr[i] = (b - sum) / matrix[i][i];
            }

            residual = normaForJakobi(xOld, xCurr);
            copy(xCurr.begin(), xCurr.end(), xOld.begin());
            cntIter++;
        }

        ans = xOld;
        return true;
    }
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

std::pair<double, double> calcNormas(std::vector<double>& ansJakobi, std::vector<double>& ansLapacke){
    double oneNorma = 0.0;
    double twoNorma = 0.0;
    double oneNormaLapacke = 0.0;
    double twoNormaLapacke = 0.0;

    for (int i = 0; i < N; ++i){
        oneNorma += fabs(ansLapacke[i] - ansJakobi[i]);
        twoNorma += (ansLapacke[i] - ansJakobi[i]) * (ansLapacke[i] - ansJakobi[i]);
        oneNormaLapacke += fabs(ansLapacke[i]);
        twoNormaLapacke += ansLapacke[i] * ansLapacke[i];
    }

    return {oneNorma / oneNormaLapacke, sqrt(twoNorma) / sqrt(twoNormaLapacke)};
}

void compareJakobiLapacke(std::vector<double>& matrixForLapacke, std::vector<double>& rightPart){
    std::vector<double> ansJakobi;
    int cntIterJacobi = 0;
    if (jakobi(ansJakobi, cntIterJacobi)){
        lapack(matrixForLapacke, rightPart);

        std::pair<double, double> normas = calcNormas(ansJakobi, rightPart);
        std::cout << "Количество итераций для нахождения решения = " << cntIterJacobi << '\n';
        std::cout << "Относительная погрешность метода Якоби по 1-норме = " << normas.first << '\n';
        std::cout << "Относительная погрешность метода Якоби по 2-норме = " << normas.second << '\n';
        std::cout << "За точное взято решение с помощью Lapacke\n\n";
    }

}

int main(){
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
   
    compareJakobiLapacke(matrixForLapacke, rightPart);

    return 0;
}
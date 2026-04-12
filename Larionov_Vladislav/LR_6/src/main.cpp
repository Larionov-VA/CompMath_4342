#include <iostream>
#include <random>

#define MAXIMUM_ELEMENT_VALUE 10
#ifndef EPSILON
#define EPSILON 1e-8
#endif


void Jacobi(int N, double** A, double* F, double* X) {
    double* TempX = new double[N]();
    double* X_prev = new double[N];
    for (int i = 0; i < N; i++) {
        X[i] = 0;
    }
    double diff;
    int iterations = 0;
    const long MAX_ITER = 10000;
    do {
        for (int i = 0; i < N; i++) {
            X_prev[i] = X[i];
        }
        for (int i = 0; i < N; i++) {
            double sum = F[i];
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    sum -= A[i][j] * X_prev[j];
                }
            }
            TempX[i] = sum / A[i][i];
        }
        for (int i = 0; i < N; i++) {
            diff = fabs(TempX[i] - X_prev[i]);
        }
        for (int i = 0; i < N; i++) {
            X[i] = TempX[i];
        }
        iterations++;
        if (iterations >= MAX_ITER) {
            std::cerr << "Превышено максимальное число итераций!\n";
            break;
        }
    } while (diff > EPSILON);
    std::cout << "Сошлось за " << iterations
        << " итераций, |x_prev - x_current| = " << diff << std::endl;
    delete[] TempX;
    delete[] X_prev;
}


double getRandomNumber() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dist(
        -MAXIMUM_ELEMENT_VALUE, MAXIMUM_ELEMENT_VALUE
    );
    return dist(gen);
}


double** createMatrix(int N) {
    double** A = new double*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
        double sum = 0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                A[i][j] = getRandomNumber();
                sum += fabs(A[i][j]);
            }
        }
        A[i][i] = sum + 10000000;
    }
    return A;
}


std::pair<double*, double*> getCorrectResponseVectors(double** A, int N) {
    double* X = new double[N];
    double* B = new double[N];
    for (int i = 0; i < N; ++i) {
        X[i] = getRandomNumber();
    }
    for (int i = 0; i < N; ++i) {
        B[i] = 0;
        for (int j = 0; j < N; ++j) {
            B[i] += A[i][j] * X[j];
        }
    }
    return {X, B};
}


void deleteMatrix(double** A, int N) {
    for (int i = 0; i < N; ++i) {
        delete[] A[i];
    }
    delete[] A;
}


int main() {
    int N;
    std::cin >> N;
    auto A = createMatrix(N);
    auto vectors = getCorrectResponseVectors(A, N);
    double* JacobiAnswer = new double[N]();
    Jacobi(N, A, vectors.second, JacobiAnswer);
    std::cout << "Решение методом Якоби:\n";
    for (int i = 0; i < N; i++) {
        std::cout << "X[" << i << "] = " << JacobiAnswer[i]
                  << " (должно быть " << vectors.first[i] << ")\n";
    }
    delete[] vectors.first;
    delete[] vectors.second;
    delete[] JacobiAnswer;
    deleteMatrix(A, N);
    return EXIT_SUCCESS;
}
#include <iostream>
#include <cmath>
#include <cstdio>
#include <random>

#define MAXIMUM_ELEMENT_VALUE 10

int Jacobi(int N, double** A, double* F, double* X, double eps) {
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
            std::cerr << "Error: iteration limit!\n";
            break;
        }
    } while (diff > eps);
    delete[] TempX;
    delete[] X_prev;

    return iterations;
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


int main(int argc, char *argv[]) {
    // iteration over precision
    FILE *precs = fopen(argv[1], "w");
    for (int i=1; i<20; i++) {
        auto A = createMatrix(25);
        auto vectors = getCorrectResponseVectors(A, 25);
        double* JacobiAnswer = new double[25]();
        int iter = Jacobi(25, A, vectors.second, JacobiAnswer, std::pow(10, -i));
        fprintf(precs, "%d %d\n", i, iter);
        delete[] vectors.first;
        delete[] vectors.second;
        delete[] JacobiAnswer;
        deleteMatrix(A, 25);
    }
    // iterations over size
    FILE *comp = fopen(argv[2], "w");
    for (int n=1; n<=1000; n++) {
        auto A = createMatrix(n);
        auto vectors = getCorrectResponseVectors(A, n);
        double* JacobiAnswer = new double[n]();
        int iter = Jacobi(n, A, vectors.second, JacobiAnswer, std::pow(10, -8));
        fprintf(comp, "%d %d\n", n, iter);
        delete[] vectors.first;
        delete[] vectors.second;
        delete[] JacobiAnswer;
        deleteMatrix(A, n);
    }
    return EXIT_SUCCESS;
}

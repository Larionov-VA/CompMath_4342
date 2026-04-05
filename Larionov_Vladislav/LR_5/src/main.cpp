#include "../utils/matrixlib.hpp"
#include <iomanip>

using type = long double;

template <class T>
std::pair<std::vector<T>, std::vector<T>>
getCorrectResponseVectors(squareMatrix<T>& M) {
    int N = M.getN();
    std::vector<T> X(N, 0), B(N, 0);
    if (!M.isSingular()) {
        for (int i = 0; i < N; ++i) {
            X[i] = M.getRandomNumber();
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                B[i] += (M.getElement(i, j) * X[j]);
            }
        }
    }
    return {X, B};
}

template <class T>
std::vector<T> GaussianElimination(squareMatrix<T>& M, std::vector<T>& B) {
    int N = M.getN();
    std::vector<T> X(N, 0);
    for (int i = 0; i < N - 1; ++i) {
        auto coords = M.getCoordOfMaxElement(i);
        if (coords.first != i) {
            M.swapRows(coords.first, i);
            std::swap(B[coords.first], B[i]);
        }
        for (int j = i + 1; j < N; ++j) {
            T factor = M.getElement(j, i) / M.getElement(i, i);
            for (int k = i; k < N; ++k) {
                T val = M.getElement(j, k) - factor * M.getElement(i, k);
                M.setElement(j, k, val);
            }
            B[j] -= factor * B[i];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        T sum = 0;
        for (int j = i + 1; j < N; ++j) {
            sum += M.getElement(i, j) * X[j];
        }
        X[i] = (B[i] - sum) / M.getElement(i, i);
    }
    return X;
}

int main() {
    int N = 0;
    std::cin >> N;
    squareMatrix<type> M = squareMatrix<type>(N);
    std::cout << "Число обусловленности матрицы: " << M.conditionNumber() << '\n';
    auto vectors = getCorrectResponseVectors<type>(M);
    std::fstream logFile("log.txt", std::ios::out);
    logFile << M << '\n';
    std::vector<type> answer = GaussianElimination<type>(M, vectors.second);
    long double maxDelta = std::numeric_limits<long double>::min();
    for (int i = 0; i < N; ++i) {
        if (answer[i] != vectors.first[i]) {
            long double currentDelta = std::fabs(vectors.first[i] - answer[i]);
            if (currentDelta > maxDelta) {
                maxDelta = currentDelta;
            }
            logFile << i << ' ' << std::fixed << std::setprecision(64) << answer[i]
            << " ~ " << vectors.first[i] << "\ndelta = " << currentDelta << '\n';
        }
        else {
            logFile << i << ' ' << std::fixed << std::setprecision(64) << answer[i] << '\n';
        }
    }
    std::cout << "Максимальная погрешность: " << maxDelta << '\n';
    logFile << M << '\n';
    return 0;
}
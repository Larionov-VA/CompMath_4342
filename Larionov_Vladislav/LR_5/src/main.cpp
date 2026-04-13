#include "Eigen/Dense"
#include "../utils/matrixlib.hpp"
#include <chrono>

#define REPEAT_COUNT 10

using type = float;

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
        int RowOfMaxElementInColumn = M.getRowOfMaxElement(i, i);
        if (RowOfMaxElementInColumn != i && RowOfMaxElementInColumn != -1) {
            M.swapRows(RowOfMaxElementInColumn, i);
            std::swap(B[RowOfMaxElementInColumn], B[i]);
        }
        auto& row_i = M.getRow(i);
        for (int j = i + 1; j < N; ++j) {
            auto& row_j = M.getRow(j);
            T factor = row_j[i] / row_i[i];
            for (int k = i; k < N; ++k)
                row_j[k] -= factor * row_i[k];
            B[j] -= factor * B[i];
        }
    }
    for (int i = N - 1; i >= 0; --i) {
        const auto& row_i = M.getRow(i);
        T sum = 0;
        for (int j = i + 1; j < N; ++j) {
            sum += row_i[j] * X[j];
        }
        X[i] = (B[i] - sum) / row_i[i];
    }
    return X;
}

int getData() {
    std::fstream log("log.txt", std::ios::out);
    std::vector<std::vector<long double>> results;
    results.reserve(200);

    for (int i = 2; i < 2000; i += 20) {
        int count = 0;
        std::vector<long double> row(5);
        row[0] = i;
        while(count++ != REPEAT_COUNT) {
            squareMatrix<type> M_original(i);
            auto vectors = getCorrectResponseVectors<type>(M_original);

            squareMatrix<type> M_for_gaussian = M_original;

            Eigen::MatrixXf eigenM(i, i);
            Eigen::VectorXf eigenB(i);
            for (int k = 0; k < i; ++k) {
                for (int j = 0; j < i; ++j) {
                    eigenM(k, j) = M_original.getElement(k, j);
                }
                eigenB(k) = vectors.second[k];
            }

            auto startEigen = std::chrono::high_resolution_clock::now();
            Eigen::VectorXf eigenX = eigenM.partialPivLu().solve(eigenB);
            auto endEigen = std::chrono::high_resolution_clock::now();
            row[1] += std::chrono::duration<double, std::milli>(endEigen - startEigen).count();

            long double err_eigen = 0.0;
            int valid_count_eigen = 0;
            for (int k = 0; k < i; ++k) {
                if (std::fabs(vectors.first[k]) > 1e-10) {
                    double d = std::fabs(vectors.first[k] - eigenX(k));
                    err_eigen += d / std::fabs(vectors.first[k]);
                    valid_count_eigen++;
                }
            }
            row[2] += (valid_count_eigen > 0) ? err_eigen / valid_count_eigen : 0.0;

            auto start = std::chrono::high_resolution_clock::now();
            std::vector<type> answer = GaussianElimination<type>(M_for_gaussian, vectors.second);
            auto end = std::chrono::high_resolution_clock::now();
            row[3] += std::chrono::duration<double, std::milli>(end - start).count();

            long double err_gaussian = 0.0;
            int valid_count_gaussian = 0;
            for (int k = 0; k < i; ++k) {
                if (std::fabs(vectors.first[k]) > 1e-10) {
                    double d = std::fabs(vectors.first[k] - answer[k]);
                    err_gaussian += d / std::fabs(vectors.first[k]);
                    valid_count_gaussian++;
                }
            }
            row[4] += (valid_count_gaussian > 0) ? err_gaussian / valid_count_gaussian : 0.0;
        }
        row[1] /= REPEAT_COUNT;
        row[2] /= REPEAT_COUNT;
        row[3] /= REPEAT_COUNT;
        row[4] /= REPEAT_COUNT;
        results.push_back(row);
    }
    for (const auto& row : results) {
        for (auto val : row) log << val << '\t';
        log << '\n';
    }
    return EXIT_SUCCESS;
}

int TestGaussian() {
    int N = 0;
    std::cin >> N;
    squareMatrix<type> M(N);
    auto vectors = getCorrectResponseVectors<type>(M);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<type> answer = GaussianElimination<type>(M, vectors.second);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Время метода GaussianElimination: " << duration.count() << " мс" << std::endl;
    long double errorNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        errorNorm += (vectors.first[i] - answer[i])*(vectors.first[i] - answer[i]);
    }
    std::cout << "||точный_X - GaussianElimination|| = " << std::sqrt(errorNorm) << std::endl;
    return EXIT_SUCCESS;
}

int CompareGaussianAndLU() {
    int N = 0;
    std::cin >> N;
    squareMatrix<type> M(N);
    auto vectors = getCorrectResponseVectors<type>(M);
    Eigen::MatrixXf eigenM(N, N);
    Eigen::VectorXf eigenB(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            eigenM(i, j) = M.getElement(i, j);
        }
        eigenB(i) = vectors.second[i];
    }

    auto startEigen = std::chrono::high_resolution_clock::now();
    Eigen::VectorXf eigenX = eigenM.partialPivLu().solve(eigenB);
    auto endEigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> durationEigen = endEigen - startEigen;
    std::cout << "Время решения Eigen: " << durationEigen.count() << " мс" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<type> answer = GaussianElimination<type>(M, vectors.second);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << "Время метода GaussianElimination: " << duration.count() << " мс" << std::endl;

    long double diffNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        diffNorm += (vectors.first[i] - eigenX(i))*(vectors.first[i] - eigenX(i));
    }
    std::cout << "||точный_X - eigen_решение|| = " << std::sqrt(diffNorm) << std::endl;

    long double errorNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        errorNorm += (vectors.first[i] - answer[i])*(vectors.first[i] - answer[i]);
    }
    std::cout << "||точный_X - GaussianElimination|| = " << std::sqrt(errorNorm) << std::endl;
    return EXIT_SUCCESS;
}

int main() {
    std::cout << "Выберете из списка:\n";
    std::cout << "1 - Получить и вывести данные для построения графика сравнения времени выполнения и точности для методов\n";
    std::cout << "2 - Сравнить метод GaussianElimination с методом LU из библиотеки Eigen для матрицы N на N\n";
    std::cout << "3 - Выполнить GaussianElimination для матрицы N на N\n";
    int choice;
    std::cin >> choice;
    switch (choice)
    {
    case 1:
        return getData();
    case 2:
        return CompareGaussianAndLU();
    case 3:
        return TestGaussian();
    }
}
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
#include "matlib.hpp"

#define REPEAT_COUNT 10
using type = double;

template <class T>
std::pair<std::vector<T>, std::vector<T>> getCorrectResponseVectors(squareMatrix<T>& M) {
    int N = M.getN();
    std::vector<T> X(N, 0), B(N, 0);
    for (int i = 0; i < N; ++i)
        X[i] = M.getRandomNumber();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            B[i] += M.getElement(i, j) * X[j];
    return {X, B};
}

template <class T>
std::vector<T> GaussianElimination(squareMatrix<T>& M, std::vector<T>& B) {
    int N = M.getN();
    std::vector<T> X(N, 0);

    for (int i = 0; i < N - 1; ++i) {
        int pivotRow = M.getRowOfMaxElement(i, i);
        if (pivotRow == -1 || std::fabs(M.getElement(pivotRow, i)) < EPSILON) {
            throw std::runtime_error("Матрица вырождена!");
        }
        if (pivotRow != i) {
            M.swapRows(pivotRow, i);
            std::swap(B[pivotRow], B[i]);
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

    if (std::fabs(M.getElement(N-1, N-1)) < EPSILON)
        throw std::runtime_error("Матрица вырождена!");

    for (int i = N - 1; i >= 0; --i) {
        const auto& row_i = M.getRow(i);
        T sum = 0;
        for (int j = i + 1; j < N; ++j)
            sum += row_i[j] * X[j];
        X[i] = (B[i] - sum) / row_i[i];
    }
    return X;
}

int getData() {
    std::ofstream log("log.txt");
    std::vector<std::vector<long double>> results;
    std::cout << "Сбор статистики (N = 2 .. 1000)\n";

    for (int N = 2; N <= 1000; N += 25) {
        long double t_eigen = 0, err_eigen = 0, t_gauss = 0, err_gauss = 0;
        int validCount = 0;

        for (int rep = 0; rep < REPEAT_COUNT; ++rep) {
            squareMatrix<type> M_orig(N);
            auto [X_exact, B] = getCorrectResponseVectors(M_orig);

            Eigen::MatrixXd eigenM(N, N);
            Eigen::VectorXd eigenB(N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j)
                    eigenM(i, j) = M_orig.getElement(i, j);
                eigenB(i) = B[i];
            }

            auto startE = std::chrono::high_resolution_clock::now();
            Eigen::VectorXd eigenX = eigenM.partialPivLu().solve(eigenB);
            auto endE = std::chrono::high_resolution_clock::now();
            t_eigen += std::chrono::duration<double, std::milli>(endE - startE).count();

            long double err_e = 0.0;
            int cnt_e = 0;
            for (int i = 0; i < N; ++i) {
                if (std::fabs(X_exact[i]) > 1e-10) {
                    err_e += std::fabs(X_exact[i] - eigenX(i)) / std::fabs(X_exact[i]);
                    ++cnt_e;
                }
            }
            err_eigen += (cnt_e > 0) ? err_e / cnt_e : 0.0;

            squareMatrix<type> M_gauss = M_orig;
            auto startG = std::chrono::high_resolution_clock::now();
            try {
                std::vector<type> X_gauss = GaussianElimination(M_gauss, B);
                auto endG = std::chrono::high_resolution_clock::now();
                t_gauss += std::chrono::duration<double, std::milli>(endG - startG).count();

                long double err_g = 0.0;
                int cnt_g = 0;
                for (int i = 0; i < N; ++i) {
                    if (std::fabs(X_exact[i]) > 1e-10) {
                        err_g += std::fabs(X_exact[i] - X_gauss[i]) / std::fabs(X_exact[i]);
                        ++cnt_g;
                    }
                }
                err_gauss += (cnt_g > 0) ? err_g / cnt_g : 0.0;
                ++validCount;
            } catch (...) {
            }
        }
        if (validCount > 0) {
            t_eigen /= REPEAT_COUNT;
            err_eigen /= REPEAT_COUNT;
            t_gauss /= REPEAT_COUNT;
            err_gauss /= REPEAT_COUNT;
            results.push_back({(long double)N, t_eigen, err_eigen, t_gauss, err_gauss});
        }
        if (N % 200 == 2) std::cout << "Обработано N = " << N << "\n";
    }

    for (auto& row : results) {
        for (auto val : row) log << val << '\t';
        log << '\n';
    }
    log.close();

    std::ofstream py("plotter.py");
    py << R"(
import matplotlib.pyplot as plt
import sys

N, t_e, e_e, t_g, e_g = [], [], [], [], []
try:
    with open('log.txt', 'r') as f:
        for line in f:
            data = list(map(float, line.split()))
            N.append(data[0]); t_e.append(data[1]); e_e.append(data[2]); t_g.append(data[3]); e_g.append(data[4])
except Exception as e:
    print('Ошибка чтения:', e)
    sys.exit(1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
ax1.plot(N, t_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax1.plot(N, t_g, 'r-s', markersize=4, label='GaussianElimination')
ax1.set_yscale('log'); ax1.grid(True, linestyle='--'); ax1.set_title('Время выполнения'); ax1.set_xlabel('N'); ax1.set_ylabel('Время, мс'); ax1.legend()

ax2.plot(N, e_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax2.plot(N, e_g, 'r-s', markersize=4, label='GaussianElimination')
ax2.set_yscale('log'); ax2.grid(True, linestyle='--'); ax2.set_title('Относительная погрешность'); ax2.set_xlabel('N'); ax2.set_ylabel('Погрешность'); ax2.legend()

plt.tight_layout()
plt.savefig('graphs.png', dpi=300)
plt.show()
)";
    py.close();
    system("python plotter.py");
    return 0;
}

int CompareGaussianAndLU() {
    int N;
    std::cout << "Введите размер матрицы N: ";
    std::cin >> N;
    squareMatrix<type> M(N);
    auto [X_exact, B] = getCorrectResponseVectors(M);

    Eigen::MatrixXd eigenM(N, N);
    Eigen::VectorXd eigenB(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) eigenM(i, j) = M.getElement(i, j);
        eigenB(i) = B[i];
    }

    auto startE = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd eigenX = eigenM.partialPivLu().solve(eigenB);
    auto endE = std::chrono::high_resolution_clock::now();
    double tEigen = std::chrono::duration<double, std::milli>(endE - startE).count();
    std::cout << "Eigen время: " << tEigen << " мс\n";

    try {
        auto startG = std::chrono::high_resolution_clock::now();
        std::vector<type> X_gauss = GaussianElimination(M, B);
        auto endG = std::chrono::high_resolution_clock::now();
        double tGauss = std::chrono::duration<double, std::milli>(endG - startG).count();
        std::cout << "Метод Гаусса время: " << tGauss << " мс\n";

        double errEigen = 0.0, errGauss = 0.0;
        for (int i = 0; i < N; ++i) {
            errEigen += (X_exact[i] - eigenX(i)) * (X_exact[i] - eigenX(i));
            errGauss += (X_exact[i] - X_gauss[i]) * (X_exact[i] - X_gauss[i]);
        }
        std::cout << "||Точный - Eigen||_2 = " << std::sqrt(errEigen) << "\n";
        std::cout << "||Точный - Gauss||_2 = " << std::sqrt(errGauss) << "\n";
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << "\n";
    }
    return 0;
}

int TestLargeMatrices() {
    std::vector<int> sizes = {1000, 10000};
    for (int N : sizes) {
        std::cout << "\n--- Тест N = " << N << " ---\n";
        try {
            squareMatrix<type> M(N);
            auto [X_exact, B] = getCorrectResponseVectors(M);
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<type> X_gauss = GaussianElimination(M, B);
            auto end = std::chrono::high_resolution_clock::now();
            double sec = std::chrono::duration<double>(end - start).count();
            std::cout << "Успех! Время: " << sec << " сек\n";
        } catch (const std::bad_alloc&) {
            std::cout << "Недостаточно памяти для N = " << N << "\n";
        } catch (const std::exception& e) {
            std::cout << "Ошибка: " << e.what() << "\n";
        }
    }
    return 0;
}

int TestHilbertMatrix() {
    int N;
    std::cout << "Введите размер матрицы Гильберта N: ";
    std::cin >> N;
    squareMatrix<type> M(N, true);
    auto [X_exact, B] = getCorrectResponseVectors(M);
    try {
        std::vector<type> X_gauss = GaussianElimination(M, B);
        double err = 0.0;
        for (int i = 0; i < N; ++i)
            err += (X_exact[i] - X_gauss[i]) * (X_exact[i] - X_gauss[i]);
        std::cout << "Евклидова норма погрешности: " << std::sqrt(err) << "\n";
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << "\n";
    }
    return 0;
}

int main() {
    std::cout << "Ввод для выполнения\n";
    std::cout << "1 - Собрать данные для графиков (N=2..1000)\n";
    std::cout << "2 - Сравнить с Eigen на одной матрице\n";
    std::cout << "3 - Стресс-тест (1000 и 10000)\n";
    std::cout << "4 - Тест на матрице Гильберта\n";
    int choice;
    std::cin >> choice;
    switch (choice) {
        case 1: return getData();
        case 2: return CompareGaussianAndLU();
        case 3: return TestLargeMatrices();
        case 4: return TestHilbertMatrix();
        default: std::cout << "Неверный ввод\n"; return 0;
    }
}
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
std::pair<std::vector<T>, std::vector<T>>
getCorrectResponseVectors(squareMatrix<T>& M) {
    int N = M.getN();
    std::vector<T> X(N, 0), B(N, 0);
    for (int i = 0; i < N; ++i) X[i] = M.getRandomNumber();
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            B[i] += M.getElement(i, j) * X[j];
    return {X, B};
}

template <class T>
std::vector<T> JacobiMethod(squareMatrix<T>& M, const std::vector<T>& B) {
    int N = M.getN();
    std::vector<T> X(N, 0.0);
    std::vector<T> X_new(N, 0.0);
    T diff;
    int iterations = 0;
    const int MAX_ITER = 10000;

    for (int i = 0; i < N; ++i)
        if (std::fabs(M.getElement(i, i)) < EPSILON)
            throw std::runtime_error("Нулевой элемент на диагонали");

    do {
        diff = 0.0;
        for (int i = 0; i < N; ++i) {
            T sum = B[i];
            for (int j = 0; j < N; ++j)
                if (i != j) sum -= M.getElement(i, j) * X[j];
            X_new[i] = sum / M.getElement(i, i);
            T delta = std::fabs(X_new[i] - X[i]);
            if (delta > diff) diff = delta;
        }
        X = X_new;
        iterations++;
        if (iterations >= MAX_ITER)
            throw std::runtime_error("Метод Якоби расходится (нет диагонального преобладания)");
    } while (diff > EPSILON);

    return X;
}

int getData() {
    std::ofstream log("log.txt");
    std::vector<std::vector<long double>> results;
    std::cout << "Сбор статистики (N от 2 до 2000, шаг 20)...\n";

    for (int N = 2; N <= 2000; N += 20) {
        long double t_eigen = 0, err_eigen = 0, t_jacobi = 0, err_jacobi = 0;
        int valid = 0;

        for (int rep = 0; rep < REPEAT_COUNT; ++rep) {
            squareMatrix<type> M_orig(N);
            M_orig.makeDiagonallyDominant();
            auto [X_exact, B] = getCorrectResponseVectors(M_orig);

            // Eigen (прямой метод)
            Eigen::MatrixXd eigenM(N, N);
            Eigen::VectorXd eigenB(N);
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) eigenM(i, j) = M_orig.getElement(i, j);
                eigenB(i) = B[i];
            }

            auto startE = std::chrono::high_resolution_clock::now();
            Eigen::VectorXd eigenX = eigenM.partialPivLu().solve(eigenB);
            auto endE = std::chrono::high_resolution_clock::now();
            t_eigen += std::chrono::duration<double, std::milli>(endE - startE).count();

            long double err_e = 0.0; int cnt_e = 0;
            for (int i = 0; i < N; ++i) {
                if (std::fabs(X_exact[i]) > 1e-12) {
                    err_e += std::fabs(X_exact[i] - eigenX(i)) / std::fabs(X_exact[i]);
                    cnt_e++;
                }
            }
            err_eigen += (cnt_e > 0) ? err_e / cnt_e : 0.0;

            squareMatrix<type> M_jac = M_orig;
            auto startJ = std::chrono::high_resolution_clock::now();
            try {
                std::vector<type> X_jac = JacobiMethod(M_jac, B);
                auto endJ = std::chrono::high_resolution_clock::now();
                t_jacobi += std::chrono::duration<double, std::milli>(endJ - startJ).count();

                long double err_j = 0.0; int cnt_j = 0;
                for (int i = 0; i < N; ++i) {
                    if (std::fabs(X_exact[i]) > 1e-12) {
                        err_j += std::fabs(X_exact[i] - X_jac[i]) / std::fabs(X_exact[i]);
                        cnt_j++;
                    }
                }
                err_jacobi += (cnt_j > 0) ? err_j / cnt_j : 0.0;
                ++valid;
            } catch (...) {
            }
        }

        if (valid > 0) {
            t_eigen /= REPEAT_COUNT;
            err_eigen /= REPEAT_COUNT;
            t_jacobi /= REPEAT_COUNT;
            err_jacobi /= REPEAT_COUNT;
            results.push_back({(long double)N, t_eigen, err_eigen, t_jacobi, err_jacobi});
        }
        if (N % 200 == 0) std::cout << "Обработано N = " << N << "\n";
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

N, t_e, e_e, t_j, e_j = [], [], [], [], []
try:
    with open('log.txt', 'r') as f:
        for line in f:
            data = list(map(float, line.split()))
            N.append(data[0]); t_e.append(data[1]); e_e.append(data[2]); t_j.append(data[3]); e_j.append(data[4])
except Exception as e:
    print('Ошибка чтения:', e)
    sys.exit(1)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,6))
ax1.plot(N, t_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax1.plot(N, t_j, 'g-s', markersize=4, label='Jacobi Method')
ax1.set_yscale('log')
ax1.grid(True, linestyle='--')
ax1.set_title('Время выполнения')
ax1.set_xlabel('Размер матрицы N')
ax1.set_ylabel('Время, мс')
ax1.legend()

ax2.plot(N, e_e, 'b-o', markersize=4, label='Eigen (partialPivLu)')
ax2.plot(N, e_j, 'g-s', markersize=4, label='Jacobi Method')
ax2.set_yscale('log')
ax2.grid(True, linestyle='--')
ax2.set_title('Относительная погрешность')
ax2.set_xlabel('Размер матрицы N')
ax2.set_ylabel('погрешность')
ax2.legend()

plt.tight_layout()
plt.savefig('graphs_jacobi.png', dpi=300)
plt.show()
)";
    py.close();
    system("python3 plotter.py");
    return 0;
}

int CompareJacobiAndLU() {
    int N;
    std::cout << "Введите размер матрицы N: ";
    std::cin >> N;
    squareMatrix<type> M(N);
    M.makeDiagonallyDominant();
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
        auto startJ = std::chrono::high_resolution_clock::now();
        std::vector<type> X_jac = JacobiMethod(M, B);
        auto endJ = std::chrono::high_resolution_clock::now();
        double tJac = std::chrono::duration<double, std::milli>(endJ - startJ).count();
        std::cout << "Метод Якоби время: " << tJac << " мс\n";

        double errEigen = 0.0, errJac = 0.0;
        for (int i = 0; i < N; ++i) {
            errEigen += (X_exact[i] - eigenX(i)) * (X_exact[i] - eigenX(i));
            errJac   += (X_exact[i] - X_jac[i])   * (X_exact[i] - X_jac[i]);
        }
        std::cout << "||Точный - Eigen||_2  = " << std::sqrt(errEigen) << "\n";
        std::cout << "||Точный - Jacobi||_2 = " << std::sqrt(errJac) << "\n";
    } catch (const std::exception& e) {
        std::cout << "Ошибка: " << e.what() << "\n";
    }
    return 0;
}

int TestLargeMatrices() {
    for (int N : {1000, 10000}) {
        std::cout << "\n--- Тест N = " << N << " ---\n";
        try {
            squareMatrix<type> M(N);
            M.makeDiagonallyDominant();
            auto [X_exact, B] = getCorrectResponseVectors(M);
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<type> X_jac = JacobiMethod(M, B);
            auto end = std::chrono::high_resolution_clock::now();
            double sec = std::chrono::duration<double>(end - start).count();
            std::cout << "Успех! Время: " << sec << " сек\n";
        } catch (const std::bad_alloc&) {
            std::cout << "Недостаточно памяти\n";
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
        std::vector<type> X_jac = JacobiMethod(M, B);
        std::cout << "Ошибка! Неожиданно сошлось!\n";
    } catch (const std::exception& e) {
        std::cout << "\n[ОЖИДАЕМАЯ ОШИБКА] " << e.what() << "\n";
        std::cout << "Матрица Гильберта не имеет диагонального преобладания → метод расходится.\n";
    }
    return 0;
}

int main() {
    std::cout << "Ввод для выполнения\n";
    std::cout << "1 - Собрать данные для графиков (N=2..2000)\n";
    std::cout << "2 - Сравнить метод Якоби с Eigen\n";
    std::cout << "3 - Стресс-тест (1000, 10000)\n";
    std::cout << "4 - Тест на матрице Гильберта (расходимость)\n";
    int choice;
    std::cin >> choice;
    switch (choice) {
        case 1: return getData();
        case 2: return CompareJacobiAndLU();
        case 3: return TestLargeMatrices();
        case 4: return TestHilbertMatrix();
        default: std::cout << "Неверный выбор\n"; return 0;
    }
}
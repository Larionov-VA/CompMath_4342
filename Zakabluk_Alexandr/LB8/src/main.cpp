#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <windows.h>

#include "Eigen/Dense"
#include "matrixlib.hpp"

#define REPEAT_COUNT 10

using type = double;

template <class T>
std::pair<std::vector<T>, std::vector<T>>
getCorrectResponseVectors(squareMatrix<T>& M) {
    int N = M.getN();
    std::vector<T> X(N, 0), B(N, 0);
    for (int i = 0; i < N; ++i) {
        X[i] = M.getRandomNumber();
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            B[i] += (M.getElement(i, j) * X[j]);
        }
    }
    return {X, B};
}

// РЕАЛИЗАЦИЯ МЕТОДА ЯКОБИ
template <class T>
std::vector<T> JacobiMethod(squareMatrix<T>& M, const std::vector<T>& B) {
    int N = M.getN();
    std::vector<T> X(N, 0.0); // Начальное приближение: нули
    std::vector<T> X_new(N, 0.0);
    T diff;
    int iterations = 0;
    const int MAX_ITER = 10000; // Защита от бесконечного цикла

    // Проверка диагонали на нули
    for (int i = 0; i < N; ++i) {
        if (std::fabs(M.getElement(i, i)) < EPSILON) {
            throw std::runtime_error("Нулевой элемент на главной диагонали! Метод неприменим.");
        }
    }

    do {
        diff = 0.0;
        for (int i = 0; i < N; ++i) {
            T sum = B[i];
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    sum -= M.getElement(i, j) * X[j];
                }
            }
            X_new[i] = sum / M.getElement(i, i);
            
            // Ищем максимальную разность по модулю (норма бесконечности)
            if (std::fabs(X_new[i] - X[i]) > diff) {
                diff = std::fabs(X_new[i] - X[i]);
            }
        }
        
        X = X_new;
        iterations++;
        
        if (iterations >= MAX_ITER) {
            throw std::runtime_error("Метод Якоби расходится (превышен лимит итераций). Нет диагонального преобладания.");
        }
    } while (diff > EPSILON);

    return X;
}

int getData() {
    std::fstream log("log.txt", std::ios::out);
    std::vector<std::vector<long double>> results;
    results.reserve(200);

    std::cout << "Начинаю сбор статистики (от N=2 до N=1000). Это займет некоторое время...\n";

    for (int i = 2; i <= 1000; i += 25) {
        int count = 0;
        std::vector<long double> row(5, 0.0);
        row[0] = i;
        while(count++ != REPEAT_COUNT) {
            squareMatrix<type> M_original(i);
            M_original.makeDiagonallyDominant(); // ВАЖНО: Делаем матрицу решаемой для Якоби
            
            auto vectors = getCorrectResponseVectors<type>(M_original);
            squareMatrix<type> M_for_jacobi = M_original;

            Eigen::MatrixXd eigenM(i, i);
            Eigen::VectorXd eigenB(i);
            for (int k = 0; k < i; ++k) {
                for (int j = 0; j < i; ++j) {
                    eigenM(k, j) = M_original.getElement(k, j);
                }
                eigenB(k) = vectors.second[k];
            }

            auto startEigen = std::chrono::high_resolution_clock::now();
            Eigen::VectorXd eigenX = eigenM.partialPivLu().solve(eigenB);
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
            try {
                std::vector<type> answer = JacobiMethod<type>(M_for_jacobi, vectors.second);
                auto end = std::chrono::high_resolution_clock::now();
                row[3] += std::chrono::duration<double, std::milli>(end - start).count();

                long double err_jacobi = 0.0;
                int valid_count_jacobi = 0;
                for (int k = 0; k < i; ++k) {
                    if (std::fabs(vectors.first[k]) > 1e-10) {
                        double d = std::fabs(vectors.first[k] - answer[k]);
                        err_jacobi += d / std::fabs(vectors.first[k]);
                        valid_count_jacobi++;
                    }
                }
                row[4] += (valid_count_jacobi > 0) ? err_jacobi / valid_count_jacobi : 0.0;
            } catch (...) {
                count--; 
                continue;
            }
        }
        row[1] /= REPEAT_COUNT;
        row[2] /= REPEAT_COUNT;
        row[3] /= REPEAT_COUNT;
        row[4] /= REPEAT_COUNT;
        results.push_back(row);
        
        if (i % 200 == 2 || i == 2) std::cout << "Обработано N = " << i << "...\n";
    }

    for (const auto& row : results) {
        for (auto val : row) log << val << '\t';
        log << '\n';
    }
    log.close();
    std::cout << "Данные сохранены в log.txt. Генерирую графики...\n";

    // Python скрипт (исправлены подписи на Якоби)
    std::fstream pyScript("plotter.py", std::ios::out);
    pyScript << "import matplotlib.pyplot as plt\nimport sys\n\n";
    pyScript << "N, t_eigen, err_eigen, t_jacobi, err_jacobi = [], [], [], [], []\n";
    pyScript << "try:\n    with open('log.txt', 'r') as f:\n        for line in f:\n            data = list(map(float, line.split()))\n";
    pyScript << "            N.append(data[0])\n            t_eigen.append(data[1])\n            err_eigen.append(data[2])\n";
    pyScript << "            t_jacobi.append(data[3])\n            err_jacobi.append(data[4])\n";
    pyScript << "except Exception as e:\n    print('Ошибка чтения:', e)\n    sys.exit(1)\n\n";
    pyScript << "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))\n\n";
    
    pyScript << "ax1.plot(N, t_eigen, 'b-o', markersize=4, label='Eigen (PartialPivLU)')\n";
    pyScript << "ax1.plot(N, t_jacobi, 'g-s', markersize=4, label='Jacobi Method')\n";
    pyScript << "ax1.set_yscale('log')\nax1.grid(True, which='both', linestyle='--', linewidth=0.5)\n";
    pyScript << "ax1.set_title('Сравнение времени выполнения')\nax1.set_xlabel('Размер матрицы N')\nax1.set_ylabel('Время, мс')\nax1.legend()\n\n";

    pyScript << "ax2.plot(N, err_eigen, 'b-o', markersize=4, label='Eigen (PartialPivLU)')\n";
    pyScript << "ax2.plot(N, err_jacobi, 'g-s', markersize=4, label='Jacobi Method')\n";
    pyScript << "ax2.set_yscale('log')\nax2.grid(True, which='both', linestyle='--', linewidth=0.5)\n";
    pyScript << "ax2.set_title('Сравнение относительной погрешности')\nax2.set_xlabel('Размер матрицы N')\nax2.set_ylabel('погрешность')\nax2.legend()\n\n";

    pyScript << "plt.tight_layout()\nplt.savefig('graphs_jacobi.png', dpi=300)\nplt.show()\n";
    pyScript.close();
    system("python plotter.py");

    return EXIT_SUCCESS;
}

int CompareJacobiAndLU() {
    int N = 0;
    std::cout << "Введите размер матрицы N: ";
    std::cin >> N;
    squareMatrix<type> M(N);
    M.makeDiagonallyDominant(); // Гарантируем сходимость
    auto vectors = getCorrectResponseVectors<type>(M);
    
    Eigen::MatrixXd eigenM(N, N);
    Eigen::VectorXd eigenB(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            eigenM(i, j) = M.getElement(i, j);
        }
        eigenB(i) = vectors.second[i];
    }

    auto startEigen = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd eigenX = eigenM.partialPivLu().solve(eigenB);
    auto endEigen = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> durationEigen = endEigen - startEigen;
    std::cout << "Время решения Eigen (LU): " << durationEigen.count() << " мс" << std::endl;

    try {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<type> answer = JacobiMethod<type>(M, vectors.second);
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double, std::milli> duration = end - start;
        std::cout << "Время метода Якоби: " << duration.count() << " мс" << std::endl;

        long double diffNorm = 0.0;
        long double errorNorm = 0.0;
        for (int i = 0; i < N; ++i) {
            diffNorm += (vectors.first[i] - eigenX(i))*(vectors.first[i] - eigenX(i));
            errorNorm += (vectors.first[i] - answer[i])*(vectors.first[i] - answer[i]);
        }
        std::cout << "||Точный_X - Eigen_X||_2 = " << std::sqrt(diffNorm) << std::endl;
        std::cout << "||Точный_X - Jacobi_X||_2 = " << std::sqrt(errorNorm) << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Ошибка метода Якоби: " << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}

int TestLargeMatrices() {
    std::vector<int> sizes = {1000, 10000};
    for (int N : sizes) {
        std::cout << "\n--- Тестирование на большой матрице N = " << N << " ---\n";
        try {
            squareMatrix<type> M(N);
            M.makeDiagonallyDominant();
            auto vectors = getCorrectResponseVectors<type>(M);
            
            std::cout << "Матрица готова. Запуск итераций Якоби..." << std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            std::vector<type> answer = JacobiMethod<type>(M, vectors.second);
            auto end = std::chrono::high_resolution_clock::now();
            
            std::chrono::duration<double> duration = end - start;
            std::cout << "УСПЕХ! Время решения: " << duration.count() << " секунд" << std::endl;
        } catch (const std::bad_alloc& e) {
            std::cout << "Ошибка: Недостаточно ОЗУ для матрицы " << N << "x" << N << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Ошибка: " << e.what() << std::endl;
        }
    }
    return EXIT_SUCCESS;
}

int TestHilbertMatrix() {
    int N = 0;
    std::cout << "Введите размер матрицы Гильберта N: ";
    std::cin >> N;
    
    squareMatrix<type> M(N, true); 
    // ВНИМАНИЕ: Мы НЕ вызываем makeDiagonallyDominant(), так как хотим показать, что Якоби упадет!
    auto vectors = getCorrectResponseVectors<type>(M);
    
    try {
        std::cout << "Запуск метода Якоби на матрице Гильберта...\n";
        std::vector<type> answer = JacobiMethod<type>(M, vectors.second);
        std::cout << "Решено успешно!" << std::endl;
    } catch (const std::exception& e) {
        std::cout << "\n[ОЖИДАЕМАЯ ОШИБКА]: " << e.what() << std::endl;
        std::cout << "Матрица Гильберта не обладает диагональным преобладанием, метод Якоби не применим.\n";
    }
    return EXIT_SUCCESS;
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    std::cout << "МЕНЮ РЕШЕНИЯ СЛУ МЕТОДОМ ЯКОБИ:\n";
    std::cout << "1 - Получить данные для построения графиков сравнения\n";
    std::cout << "2 - Сравнить метод Якоби с библиотекой Eigen (N на N)\n";
    std::cout << "3 - СТРЕСС-ТЕСТ: матрицы большого размера (1000 и 10000)\n";
    std::cout << "4 - СПЕЦ. МАТРИЦЫ: Тестирование на матрице Гильберта (Тест на расходимость)\n";
    
    int choice;
    std::cin >> choice;
    switch (choice)
    {
    case 1: return getData();
    case 2: return CompareJacobiAndLU();
    case 3: return TestLargeMatrices();
    case 4: return TestHilbertMatrix();
    default: return 0;
    }
}
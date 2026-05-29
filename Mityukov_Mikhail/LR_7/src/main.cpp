#include "gauss.hpp"

#include <iostream>
#include <vector>

// Формат ввода:
//
// n
// a11 a12 ... a1n
// a21 a22 ... a2n
// ...
// an1 an2 ... ann
// b1 b2 ... bn

int main()
{
    try {
        int n;
        std::cin >> n;

        if (!std::cin || n <= 0) {
            std::cerr << "Ошибка: некорректный размер матрицы.\n";
            return 1;
        }

        std::vector<std::vector<double>> A(n, std::vector<double>(n));
        std::vector<double> b(n);

        // Считываем матрицу коэффициентов A.
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                std::cin >> A[i][j];
                if (!std::cin) {
                    std::cerr << "Ошибка: не удалось считать элемент матрицы A.\n";
                    return 1;
                }
            }
        }

        // Считываем вектор правых частей b.
        for (int i = 0; i < n; ++i) {
            std::cin >> b[i];
            if (!std::cin) {
                std::cerr << "Ошибка: не удалось считать элемент вектора b.\n";
                return 1;
            }
        }

        // Создаём решатель.
        GaussianSolver solver(1e-12);

        // Решаем систему.
        const GaussResult result = solver.solve(A, b);

        // Печатаем результат.
        printGaussResult(result);
    }
    catch (const std::exception& ex) {
        std::cerr << "Исключение: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
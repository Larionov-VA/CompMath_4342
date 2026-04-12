#include "gauss_solver.h"

#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    // Формат ввода:
    // m n
    // a11 a12 ... a1n
    // ...
    // am1 am2 ... amn
    // b1 b2 ... bm
    std::size_t m = 0;
    std::size_t n = 0;

    std::cout << "Введите m и n: ";
    if (!(std::cin >> m >> n)) {
        std::cerr << "Ошибка ввода: не удалось прочитать размеры матрицы.\n";
        return 1;
    }

    std::vector<std::vector<double>> A(m, std::vector<double>(n, 0.0));
    std::vector<double> b(m, 0.0);

    std::cout << "Введите матрицу A (" << m << " строк, " << n << " столбцов):\n";
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            if (!(std::cin >> A[i][j])) {
                std::cerr << "Ошибка ввода: не удалось прочитать A[" << i << "][" << j << "].\n";
                return 1;
            }
        }
    }

    std::cout << "Введите вектор b (" << m << " значений):\n";
    for (std::size_t i = 0; i < m; ++i) {
        if (!(std::cin >> b[i])) {
            std::cerr << "Ошибка ввода: не удалось прочитать b[" << i << "].\n";
            return 1;
        }
    }

    // Вызываем решатель. Для больших задач передаем данные через std::move,
    // чтобы не делать лишних копий больших структур.
    const GaussResult result = solveGaussian(std::move(A), std::move(b));

    // Печатаем текстовое описание результата, ранг и признак невырожденности.
    std::cout << result.message << "\n";
    std::cout << "Ранг: " << result.rank << "\n";
    std::cout << "Невырожденность: " << (result.nonsingular ? "да" : "нет") << "\n";

    // Для совместных систем печатаем найденный вектор (для бесконечного множества решений
    // это один из допустимых частных вариантов).
    if (result.type == SolutionType::Unique || result.type == SolutionType::Infinite) {
        std::cout << std::fixed << std::setprecision(8);
        for (std::size_t i = 0; i < result.x.size(); ++i) {
            std::cout << "x[" << i << "] = " << result.x[i] << "\n";
        }
    }

    return 0;
}

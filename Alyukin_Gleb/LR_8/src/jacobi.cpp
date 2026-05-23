// g++ -O2 -std=c++17 jacobi.cpp -o jacobi
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <algorithm>

enum class SolveStatus {
    Unique,     // Сходимость достигнута
    Infinite,   // Бесконечно много (в итерационных методах встречается редко)
    NoSolution  // Расходится или деление на 0
};

struct SparseEntry {
    int col;
    double val;
};

struct JacobiResult {
    SolveStatus status;
    std::vector<double> x;
    int iterations;
    double final_error;
    int rank;           // Добавлено для совместимости
    bool nonsingular;   // Добавлено для совместимости
};


using SparseMatrix = std::vector<std::vector<SparseEntry>>;

JacobiResult jacobiSolve(std::vector<std::vector<double>> A, std::vector<double> b, int n, double eps, int max_iter) {
    JacobiResult res;
    res.iterations = 0;
    res.x.assign(n, 0.0);
    res.rank = 0;
    res.nonsingular = false;

    // 1. Частичный выбор ведущего элемента, ранг
    for (int j = 0; j < n; ++j) {
        int max_row = j;
        double max_val = std::abs(A[j][j]);
        for (int i = j + 1; i < n; ++i) {
            if (std::abs(A[i][j]) > max_val) {
                max_val = std::abs(A[i][j]);
                max_row = i;
            }
        }

        if (max_val > 1e-15) {
            res.rank++; // Считаем ранг по количеству ненулевых диагональных элементов после выбора
            if (max_row != j) {
                std::swap(A[j], A[max_row]);
                std::swap(b[j], b[max_row]);
            }
        }
    }

    res.nonsingular = (res.rank == n);

    // Если ранг не полный, метод Якоби обычно не применим
    if (res.rank < n) {
        // Проверка: если система совместна при неполном ранге — бесконечно много
        // Но для итерационных методов это критическая ошибка (деление на 0 на диагонали)
        res.status = (res.rank < n) ? SolveStatus::Infinite : SolveStatus::NoSolution;
        // На практике Якоби требует полных данных на диагонали:
        for(int i=0; i<n; ++i) if(std::abs(A[i][i]) < 1e-15) {
            res.status = SolveStatus::NoSolution;
            return res;
        }
    }

    // 2. Перевод в разреженный формат
    SparseMatrix sparseA(n);
    std::vector<double> diag(n);
    for (int i = 0; i < n; ++i) {
        diag[i] = A[i][i];
        for (int j = 0; j < n; ++j) {
            if (std::abs(A[i][j]) > 1e-18) {
                sparseA[i].push_back({j, A[i][j]});
            }
        }
    }

    // 3. Итерационный процесс
    std::vector<double> x_new(n, 0.0);
    for (int k = 0; k < max_iter; ++k) {
        res.iterations++;
        double current_max_diff = 0.0;

        for (int i = 0; i < n; ++i) {   // Цикл по каждой переменной x_i
            double sum = b[i];
            for (const auto& entry : sparseA[i]) {  // Идем по ненулевым элементам
                if (entry.col != i) sum -= entry.val * res.x[entry.col]; // Если не на диагонали, выч A_ij * x_j
            }
            x_new[i] = sum / diag[i];   // делим на диагональный элемент A_ii
            current_max_diff = std::max(current_max_diff, std::abs(x_new[i] - res.x[i]));
        }

        res.x = x_new;
        res.final_error = current_max_diff;

        if (current_max_diff < eps) {
            res.status = SolveStatus::Unique;
            return res;
        }

        if (current_max_diff > 1e15) {
            res.status = SolveStatus::NoSolution;
            return res;
        }
    }

    res.status = SolveStatus::NoSolution; // Если не сошлось за лимит - считаем неудачей
    return res;
}

// int main() {
//     int n;
//     std::cout << "Введите размерность квадратной матрицы n: ";
//     if (!(std::cin >> n) || n <= 0) return 1;

//     std::vector<std::vector<double>> A(n, std::vector<double>(n));
//     std::vector<double> b(n);
//     // double eps = 1e-7;
//     std::cout << "Выберите способ заполнения:\n1) Случайно (с диагональным преобладанием)\n2) Вручную\nВыбор: ";
//     int choice;
//     std::cin >> choice;

//     if (choice == 1) {
//         std::random_device rd;
//         std::mt19937 gen(rd());
//         std::uniform_real_distribution<double> dist(-10.0, 10.0);
//         for (int i = 0; i < n; ++i) {
//             double row_sum = 0;
//             for (int j = 0; j < n; ++j) {
//                 A[i][j] = dist(gen);
//                 row_sum += std::abs(A[i][j]);
//             }
//             // Обеспечиваем сходимость (диагональное преобладание)
//             A[i][i] = (row_sum + std::abs(dist(gen))) * (dist(gen) > 0 ? 1 : -1);
//             b[i] = dist(gen);
//         }

//         // Вывод сгенерированной матрицы
//         std::cout << "\nСгенерированная матрица [A | b]:\n";
//         for (int i = 0; i < n; ++i) {
//             for (int j = 0; j < n; ++j) {
//                 std::cout << std::setw(10) << std::fixed << std::setprecision(3) << A[i][j];
//             }
//             std::cout << " | " << std::setw(10) << b[i] << "\n";
//         }

//     } else {
//         std::cout << "Введите элементы матрицы A и вектора b:\n";
//         for (int i = 0; i < n; ++i) {
//             for (int j = 0; j < n; ++j) std::cin >> A[i][j];
//             std::cin >> b[i];
//         }
//     }

//     JacobiResult res = jacobiSolve(A, b, n, 1e-7, 1000);

//     // --- Вывод ---
//     std::cout << "\n--- Результат ---\n";
//     switch (res.status) {
//         case SolveStatus::Unique:
//             std::cout << "статус: единственное решение\n";  break;
//         case SolveStatus::Infinite:
//             std::cout << "статус: бесконечно много решений\n";  break;
//         case SolveStatus::NoSolution:
//             std::cout << "статус: нет решений (или метод разошелся)\n";   break;
//     }

//     std::cout << "ранг: "          << res.rank                          << "\n";
//     std::cout << "невырожденная: " << (res.nonsingular ? "да" : "нет")  << "\n";
//     std::cout << "итераций: "      << res.iterations                    << "\n";
//     std::cout << "погрешность: "   << std::scientific << std::setprecision(2) << res.final_error << "\n";

//     if (res.status != SolveStatus::NoSolution) {
//         std::cout << std::fixed << std::setprecision(6) << "вектор x:";
//         for (double val : res.x) std::cout << " " << val;
//         std::cout << "\n";
//     }

//     return 0;
// }

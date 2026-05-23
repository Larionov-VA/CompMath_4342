// g++ -O2 -std=c++17 -DNO_MAIN jacobi.cpp test.cpp -o test
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <sstream>
#include <functional>
#include <numeric>

#include "jacobi.h"

struct TestResult {
    std::string name;
    bool        passed;
    std::string fail_reason;
    double      elapsed_ms;
    size_t      memory_bytes;
    int         iterations;
    double      error;
    SolveStatus status;
};

struct DualLog {
    std::ofstream file;
    DualLog(const std::string& fname) : file(fname) {}
    template<typename T>
    DualLog& operator<<(const T& v) {
        std::cout << v;
        file << v;
        return *this;
    }
} log_("jacobi_test_results.txt");

double measure_ms(std::function<void()> fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void print_result(const TestResult& r) {
    log_ << (r.passed ? "[ OK ]" : "[FAIL]") << "  "
         << std::left  << std::setw(45) << r.name
         << std::right << std::setw(10) << std::fixed << std::setprecision(2) << r.elapsed_ms << " ms"
         << " | iter: " << std::setw(5) << r.iterations
         << " | err: " << std::setw(10) << std::scientific << std::setprecision(2) << r.error
         << "\n";
}

TestResult run_test(const std::string& name, std::vector<std::vector<double>> A, std::vector<double> b,
                     double eps, int max_iter) {
    TestResult tr;
    tr.name = name;

    JacobiResult res;
    tr.elapsed_ms = measure_ms([&]{
        res = jacobiSolve(A, b, A.size(), eps, max_iter);
    });

    tr.iterations = res.iterations;
    tr.error = res.final_error;
    tr.status = res.status;

    // СТРОГАЯ ЛОГИКА: прошел только если Unique
    tr.passed = (res.status == SolveStatus::Unique);
    if (!tr.passed) tr.fail_reason = "Метод не сошелся";

    return tr;
}

// --- ТЕСТЫ ---
int main() {
    log_ << "Тестирование метода Якоби (Проверка сходимости)\n";
    log_ << "---------------------------------------------------\n";

    std::vector<TestResult> results;

    // // Успешные случаи
    // // 1. Базовый уникальный случай 2x2
    // // 2x + y = 5
    // // x + 3y = 10  => x=1, y=3
    // // (Есть диагональное преобладание: 2 > 1 и 3 > 1)
    // results.push_back(run_test("Base Unique 2x2", {{2, 1}, {1, 3}}, {5, 10}, 1e-8, 1000));
    // // Диагональное преобладание (Сходится всегда)
    // results.push_back(run_test("Диагональное преобладание", {{10, -1, 2}, {-1, 11, -1}, {2, -1, 10}}, {6, 25, -11}, 1e-8, 5000));
    // // Матрица без диагонального преобладания (но с перестановкой строк)
    // results.push_back(run_test("Pivoting (Перестановка строк)", {{1, 5}, {4, 1}}, {11, 6}, 1e-8, 5000));
    // results.push_back(run_test("Матрица 1000x1000", {}, {}, 1e-8, 5000)); // Заглушка

    // results.pop_back();

    // // Плохие случаи
    // // Квадратная 3x3 (требует перестановки)
    // //  2x +  y -  z =  8
    // // -3x -  y + 2z = -11
    // // -2x +  y + 2z = -3  => x=2, y=3, z=-1
    // results.push_back(run_test("Квадратная 3x3 (Нет преобладания)", {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}}, {8, -11, -3}, 1e-8, 500));

    // // Вырожденная 3x3 (rank=2) - Отсутствие диагонального преобладания
    // results.push_back(run_test("Вырожденная 3x3 (rank=2)", {{1, 2, 3}, {4, 5, 6}, {5, 7, 9}}, {6, 15, 21}, 1e-8, 100));
    // results.push_back(run_test("Hilbert 5x5 (Плохая обусловленность)",
    //     {{1, 0.5, 0.33, 0.25, 0.2}, {0.5, 0.33, 0.25, 0.2, 0.16}, {0.33, 0.25, 0.2, 0.16, 0.14}, {0.25, 0.2, 0.16, 0.14, 0.12}, {0.2, 0.16, 0.14, 0.12, 0.1}},
    //     {1, 1, 1, 1, 1}, 1e-8, 500));

    // // Нагрузочные тесты
    // int n1 = 1000;
    // std::vector<std::vector<double>> A1(n1, std::vector<double>(n1, 0.0));
    // std::vector<double> b1(n1, 1.0);
    // for(int i=0; i<n1; ++i) { A1[i][i] = 4.0; if(i>0) A1[i][i-1] = -1.0; if(i<n1-1) A1[i][i+1] = -1.0; }
    // results.push_back(run_test("Матрица 1000x1000", A1, b1, 1e-8, 5000));

    int n2 = 12000;
    std::vector<std::vector<double>> A2(n2, std::vector<double>(n2, 0.0));
    std::vector<double> b2(n2, 1.0);
    for(int i=0; i<n2; ++i) { A2[i][i] = 4.0; if(i>0) A2[i][i-1] = -1.0; if(i<n2-1) A2[i][i+1] = -1.0; }
    results.push_back(run_test("Матрица 12000x12000", A2, b2, 1e-3, 5000));

    // Вывод результатов
    int passed_count = 0;
    for (const auto& r : results) {
        print_result(r);
        if (r.passed) passed_count++;
    }

    return 0;
}

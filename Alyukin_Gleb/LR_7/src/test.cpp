// исключив main
// g++ -O2 -std=c++17 -DNO_MAIN gauss.cpp test.cpp -o test

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <sstream>
#include <functional>

#include "gauss.h"


struct TestResult {
    std::string name;           // название теста
    bool        passed;         // OK / FAIL
    std::string fail_reason;    // причина провала (пусто при OK)
    double      elapsed_ms;     // время работы, мс
    size_t      memory_bytes;   // оценка расхода памяти, байт
    int         rank;           // ранг, полученный решателем
    bool        nonsingular;    // признак невырожденности
    SolveStatus status;         // статус решения
};

// Глобальный лог: консоль + файл одновременно
struct DualLog {
    std::ofstream file;
    DualLog(const std::string& fname) : file(fname) {}

    template<typename T>
    DualLog& operator<<(const T& v) {
        std::cout << v;
        file     << v;
        return *this;
    }
} log_("test_results.txt");

// Замер времени — возвращает миллисекунды
template<typename Fn>
double measure_ms(Fn&& fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// Оценка памяти для матрицы m×n + вектора m (в байтах, без накладных расходов)
size_t estimate_memory(int m, int n) {
    // матрица A: m*n double + m указателей vector
    // вектор b: m double
    // pivot_row: n int
    // вектор x: n double
    return static_cast<size_t>(m) * n * sizeof(double)
         + static_cast<size_t>(m) * sizeof(std::vector<double>)  // заголовки строк
         + static_cast<size_t>(m) * sizeof(double)               // b
         + static_cast<size_t>(n) * sizeof(int)                  // pivot_row
         + static_cast<size_t>(n) * sizeof(double);              // x
}

// Невязка: max |A*x - b|
double residual(const std::vector<std::vector<double>>& A,
                const std::vector<double>& b,
                const std::vector<double>& x)
{
    int m = static_cast<int>(A.size());
    int n = static_cast<int>(x.size());
    double max_err = 0.0;
    for (int i = 0; i < m; ++i) {
        double s = 0.0;
        for (int j = 0; j < n; ++j) s += A[i][j] * x[j];
        max_err = std::max(max_err, std::abs(s - b[i]));
    }
    return max_err;
}

// Проверка x против известного эталона
double solution_error(const std::vector<double>& got,
                      const std::vector<double>& expected)
{
    double err = 0.0;
    for (size_t i = 0; i < got.size(); ++i)
        err = std::max(err, std::abs(got[i] - expected[i]));
    return err;
}

// Печать одной строки итогового отчёта
void print_result(const TestResult& r) {
    // Заголовок печатается один раз перед циклом, здесь только строка результата
    log_ << (r.passed ? "[ OK ]" : "[FAIL]") << "  "
         << std::left  << std::setw(42) << r.name
         << std::right << std::setw(10) << std::fixed << std::setprecision(3) << r.elapsed_ms << " ms"
         << "  |  mem " << std::setw(10) << r.memory_bytes
         << "  |  rank " << std::setw(3) << r.rank
         << "  |  nonsing " << (r.nonsingular ? "Y" : "N");
    log_ << "\n";
}


// ТЕСТЫ

// Вспомогательная обёртка: запускает тест, замеряет время и память
TestResult run_test(
    const std::string& name,
    std::vector<std::vector<double>> A,
    std::vector<double> b,
    // принимает результат и возвращает "" (ОК) или описание ошибки
    std::function<std::string(const GaussResult&,
                              const std::vector<std::vector<double>>&,
                              const std::vector<double>&)> check)
{
    TestResult tr;
    tr.name         = name;
    tr.memory_bytes = estimate_memory(static_cast<int>(A.size()),
                                      static_cast<int>(A[0].size()));

    GaussResult res;
    tr.elapsed_ms = measure_ms([&]{ res = gaussSolve(A, b,
                                        static_cast<int>(A.size()),
                                        static_cast<int>(A[0].size())); });
    tr.rank        = res.rank;
    tr.nonsingular = res.nonsingular;
    tr.status      = res.status;
    tr.fail_reason = check(res, A, b);
    tr.passed      = tr.fail_reason.empty();
    return tr;
}


// Тест 1. Единственное решение 2×2
// 2x + y  = 5
// x + 3y = 10   →   x=1, y=3

TestResult test_unique_2x2() {
    std::vector<std::vector<double>> A = {{2, 1}, {1, 3}};
    std::vector<double> b = {5, 10};
    std::vector<double> expected = {1.0, 3.0};

    return run_test("Unique 2x2 (x=1, y=3)", A, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique)    return "ожидался Unique";
            if (r.rank != 2)                        return "ожидался rank=2";
            if (!r.nonsingular)                     return "ожидалась невырожденность";
            double err = solution_error(r.x, expected);
            if (err > 1e-9) return "погрешность решения " + std::to_string(err);
            return "";
        });
}


// Тест 2. Единственное решение 3×3 (классический пример)
//   2x +  y -  z =  8
//  -3x -  y + 2z = -11
//  -2x +  y + 2z = -3   →   x=2, y=3, z=-1
// ─────────────────────────────────────────────────────────────────────────────
TestResult test_unique_3x3() {
    std::vector<std::vector<double>> A = {
        { 2,  1, -1},
        {-3, -1,  2},
        {-2,  1,  2}
    };
    std::vector<double> b = {8, -11, -3};
    std::vector<double> expected = {2.0, 3.0, -1.0};

    return run_test("Unique 3x3 (x=2, y=3, z=-1)", A, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique)    return "ожидался Unique";
            if (r.rank != 3)                        return "ожидался rank=3";
            if (!r.nonsingular)                     return "ожидалась невырожденность";
            double err = solution_error(r.x, expected);
            if (err > 1e-9) return "погрешность решения " + std::to_string(err);
            return "";
        });
}


// Тест 3. Нет решений
//   x + y = 1
//   x + y = 2   →   противоречие

TestResult test_no_solution() {
    std::vector<std::vector<double>> A = {{1, 1}, {1, 1}};
    std::vector<double> b = {1, 2};

    return run_test("NoSolution (противоречие)", A, b,
        [](const GaussResult& r, auto&, auto&) -> std::string {
            if (r.status != SolveStatus::NoSolution) return "ожидался NoSolution";
            if (r.rank != 1)                         return "ожидался rank=1";
            return "";
        });
}

// Тест 4. Нет решений (3 уравнения, явное противоречие в третьем)
//   x + 2y + 3z = 6
//   2x + 4y + 6z = 12
//   x + 2y + 3z = 7   →   последнее уравнение противоречит первому

TestResult test_no_solution_3x3() {
    std::vector<std::vector<double>> A = {
        {1, 2, 3},
        {2, 4, 6},
        {1, 2, 3}
    };
    std::vector<double> b = {6, 12, 7};

    return run_test("NoSolution 3x3 (строка 0≠c)", A, b,
        [](const GaussResult& r, auto&, auto&) -> std::string {
            if (r.status != SolveStatus::NoSolution) return "ожидался NoSolution";
            return "";
        });
}


// Тест 5. Бесконечно много решений (недоопределённая система)
//   x + 2y + 3z = 6
//   2x + 4y + 6z = 12   →   rank=1, z и y свободны
//
//   Проверяем, что частное решение (свободные = 0) удовлетворяет системе.

TestResult test_infinite() {
    std::vector<std::vector<double>> A = {
        {1, 2, 3},
        {2, 4, 6}
    };
    std::vector<double> b = {6, 12};

    return run_test("Infinite 2x3 (rank=1)", A, b,
        [](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Infinite) return "ожидался Infinite";
            if (r.rank != 1)                       return "ожидался rank=1";
            // Частное решение должно удовлетворять системе
            double res = residual(A2, b2, r.x);
            if (res > 1e-9) return "невязка частного решения " + std::to_string(res);
            return "";
        });
}


// Тест 6. Бесконечно много решений — ранг 2, матрица 3×4
TestResult test_infinite_3x4() {
    // Строим систему: первые два уравнения независимы, третье = сумма первых двух
    std::vector<std::vector<double>> A = {
        {1, 0, 1, 2},
        {0, 1, 2, 1},
        {1, 1, 3, 3}   // = строка0 + строка1
    };
    std::vector<double> b = {3, 4, 7};  // b[2] = b[0]+b[1]

    return run_test("Infinite 3x4 (rank=2)", A, b,
        [](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Infinite) return "ожидался Infinite";
            if (r.rank != 2)                       return "ожидался rank=2";
            double res = residual(A2, b2, r.x);
            if (res > 1e-9) return "невязка частного решения " + std::to_string(res);
            return "";
        });
}

// Тест 7. Матрица Гильберта 5×5 (плохо обусловленная)
//   H[i][j] = 1 / (i + j + 1)
//   Строим b = H * x_exact, где x_exact = (1,1,1,1,1), затем решаем.
//   Принимаем, если невязка < 1e-6 (число обусловленности ~5e5).
TestResult test_hilbert_5x5() {
    const int N = 5;
    std::vector<std::vector<double>> H(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            H[i][j] = 1.0 / (i + j + 1);

    std::vector<double> x_exact(N, 1.0);
    std::vector<double> b(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            b[i] += H[i][j] * x_exact[j];

    return run_test("Hilbert 5x5 (точность, cond~5e5)", H, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique) return "ожидался Unique";
            if (!r.nonsingular)                  return "ожидалась невырожденность";
            double err = solution_error(r.x, x_exact);
            // Матрица Гильберта 5×5 — число обусловленности ~5e5,
            // допускаем погрешность до 1e-6
            if (err > 1e-6) return "погрешность " + std::to_string(err) + " > 1e-6";
            return "";
        });
}

// Тест 8. Матрица Гильберта 8×8 (сильно плохо обусловленная)
//   Число обусловленности ~1.5e10 — проверяем лишь совместность и невязку.
TestResult test_hilbert_8x8() {
    const int N = 8;
    std::vector<std::vector<double>> H(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            H[i][j] = 1.0 / (i + j + 1);

    std::vector<double> x_exact(N, 1.0);
    std::vector<double> b(N, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            b[i] += H[i][j] * x_exact[j];

    return run_test("Hilbert 8x8 (сильный cond~1.5e10)", H, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique) return "ожидался Unique";
            // При таком числе обусловленности допускаем невязку до 1e-4
            double res = residual(A2, b2, r.x);
            if (res > 1e-4) return "невязка " + std::to_string(res) + " > 1e-4";
            return "";
        });
}

// Тест 9. Диагональная матрица 1000×1000
//   A = diag(1, 2, 3, ..., 1000)
//   b[i] = (i+1) * (i+1)   →   x[i] = i+1
//   Проверяем время, правильность, невязку.
TestResult test_diagonal_1000x1000() {
    const int N = 1000;
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> b(N), x_exact(N);
    for (int i = 0; i < N; ++i) {
        A[i][i] = static_cast<double>(i + 1);
        x_exact[i] = static_cast<double>(i + 1);
        b[i] = A[i][i] * x_exact[i];   // (i+1)^2
    }

    return run_test("Diagonal 1000x1000 (скорость)", A, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique) return "ожидался Unique";
            if (r.rank != N) return "ожидался rank=1000";
            if (!r.nonsingular) return "ожидалась невырожденность";
            double err = solution_error(r.x, x_exact);
            if (err > 1e-6) return "погрешность " + std::to_string(err);
            return "";
        });
}


// Тест 10. Нулевая матрица с ненулевым b (крайний случай)
//   Все A[i][j] = 0, b = (1, 0) → нет решений
TestResult test_zero_matrix() {
    std::vector<std::vector<double>> A = {{0, 0}, {0, 0}};
    std::vector<double> b = {1, 0};

    return run_test("Zero matrix, b≠0 (NoSolution)", A, b,
        [](const GaussResult& r, auto&, auto&) -> std::string {
            if (r.status != SolveStatus::NoSolution) return "ожидался NoSolution";
            if (r.rank != 0)                         return "ожидался rank=0";
            return "";
        });
}

// Тест 11. Нулевая матрица с нулевым b (бесконечно много решений)
//   Все A[i][j] = 0, b = (0, 0) → бесконечно много
TestResult test_zero_matrix_zero_b() {
    std::vector<std::vector<double>> A = {{0, 0}, {0, 0}};
    std::vector<double> b = {0, 0};

    return run_test("Zero matrix, b=0  (Infinite)", A, b,
        [](const GaussResult& r, auto&, auto&) -> std::string {
            if (r.status != SolveStatus::Infinite) return "ожидался Infinite";
            if (r.rank != 0)                       return "ожидался rank=0";
            return "";
        });
}

// Тест 12. Переопределённая совместная система (3 уравнения, 2 неизвестных)
//   x + y = 3
//   2x - y = 0    →   x=1, y=2
//   x - y = -1
TestResult test_overdetermined() {
    std::vector<std::vector<double>> A = {
        {1,  1},
        {2, -1},
        {1, -1}
    };
    std::vector<double> b = {3, 0, -1};
    std::vector<double> expected = {1.0, 2.0};

    return run_test("Overdetermined 3x2, Unique (x=1,y=2)", A, b,
        [&](const GaussResult& r, auto& A2, auto& b2) -> std::string {
            if (r.status != SolveStatus::Unique) return "ожидался Unique";
            if (r.rank != 2)                     return "ожидался rank=2";
            double err = solution_error(r.x, expected);
            if (err > 1e-9) return "погрешность " + std::to_string(err);
            return "";
        });
}


// main: запуск всех тестов и печать итогового отчёта
int main() {
    // Собираем все тесты в список
    std::vector<std::function<TestResult()>> tests = {
        test_unique_2x2,
        test_unique_3x3,
        test_no_solution,
        test_no_solution_3x3,
        test_infinite,
        test_infinite_3x4,
        test_hilbert_5x5,
        test_hilbert_8x8,
        test_diagonal_1000x1000,
        test_zero_matrix,
        test_zero_matrix_zero_b,
        test_overdetermined,
    };

    log_ << " Тесты метода Гаусса с частичным выбором ведущего элемента\n";

    std::vector<TestResult> results;
    results.reserve(tests.size());

    for (auto& t : tests) {
        TestResult r = t();
        results.push_back(r);
        print_result(r);
    }

    // Итоговая статистика
    int ok   = 0;
    int fail = 0;
    for (auto& r : results) r.passed ? ++ok : ++fail;

    log_ << "\nИтого: " << ok << " OK,  " << fail << " FAIL"
         << "  (из " << results.size() << " тестов)\n";

    if (fail > 0) {
        log_ << "\nПровалившиеся тесты:\n";
        for (auto& r : results)
            if (!r.passed)
                log_ << "  * " << r.name << " → " << r.fail_reason << "\n";
    }

    log_ << "Сохранено в test_results.txt\n";

    return fail == 0 ? 0 : 1;
}
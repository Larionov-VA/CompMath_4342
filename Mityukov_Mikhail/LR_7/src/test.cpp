#include "gauss.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#define RUN_VERY_LARGE_TESTS
 
// Простой мини-фреймворк для тестирования.
// Мы не используем сторонние библиотеки,
// чтобы файл можно было легко собрать "как есть".

// Счётчики тестов.
static int g_totalTests = 0;
static int g_passedTests = 0;

 
// Проверка условия.
// Если условие ложно - бросаем исключение.
 
void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

 
// Проверка чисел на приблизительное равенство.
 
bool approxEqual(double a, double b, double eps = 1e-9)
{
    return std::fabs(a - b) <= eps;
}

 
// Проверка двух векторов на приблизительное равенство.
 
bool vectorApproxEqual(const std::vector<double>& a,
                       const std::vector<double>& b,
                       double eps = 1e-9)
{
    if (a.size() != b.size()) {
        return false;
    }

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!approxEqual(a[i], b[i], eps)) {
            return false;
        }
    }

    return true;
}

 
// Умножение матрицы A на вектор x.
// Нужно для проверки невязки в тестах.
 
std::vector<double> multiply(const std::vector<std::vector<double>>& A,
                             const std::vector<double>& x)
{
    const int n = static_cast<int>(A.size());
    std::vector<double> result(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < static_cast<int>(A[i].size()); ++j) {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}

 
// Максимальная норма невязки:
// max_i |(A*x)_i - b_i|
 
double residualInfNorm(const std::vector<std::vector<double>>& A,
                       const std::vector<double>& x,
                       const std::vector<double>& b)
{
    const std::vector<double> Ax = multiply(A, x);

    double maxValue = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        maxValue = std::max(maxValue, std::fabs(Ax[i] - b[i]));
    }
    return maxValue;
}

 
// Генерация единичной матрицы размера n.
 
std::vector<std::vector<double>> makeIdentityMatrix(int n)
{
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        A[i][i] = 1.0;
    }
    return A;
}

 
// Генерация диагональной матрицы.
// Диагональ: 1, 2, 3, ..., n
//
// Такой тест полезен для больших размеров,
// потому что система хорошо обусловлена,
// а ожидаемое решение легко проверить.
 
std::vector<std::vector<double>> makeDiagonalMatrix(int n)
{
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        A[i][i] = static_cast<double>(i + 1);
    }
    return A;
}

 
// Генерация матрицы Гильберта:
// H(i, j) = 1 / (i + j + 1), если индексация с 0.
//
// Матрицы Гильберта - классический пример плохо обусловленных
// матриц. Они нужны для проверки численной устойчивости метода.
 
std::vector<std::vector<double>> makeHilbertMatrix(int n)
{
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            H[i][j] = 1.0 / static_cast<double>(i + j + 1);
        }
    }
    return H;
}

 
// Удобный запуск одного теста.
 
template <typename Func>
void runTest(const std::string& testName, Func testFunc)
{
    ++g_totalTests;

    try {
        testFunc();
        ++g_passedTests;
        std::cout << "[OK]   " << testName << "\n";
    }
    catch (const std::exception& ex) {
        std::cout << "[FAIL] " << testName << " : " << ex.what() << "\n";
    }
    catch (...) {
        std::cout << "[FAIL] " << testName << " : неизвестная ошибка\n";
    }
}

 
// Базовый тест: обычная система 2x2 с единственным решением.
 
void testBasicUnique()
{
    std::vector<std::vector<double>> A = {
        {2.0, 1.0},
        {5.0, 7.0}
    };
    std::vector<double> b = {11.0, 13.0};

    GaussianSolver solver;
    const GaussResult result = solver.solve(A, b);

    require(result.type == SolutionType::Unique, "Ожидалось единственное решение.");
    require(!result.isDegenerate, "Матрица не должна быть вырожденной.");

    const std::vector<double> expected = {64.0 / 9.0, -29.0 / 9.0};
    require(vectorApproxEqual(result.x, expected, 1e-9), "Неверное решение 2x2.");
}

 
// Тест, где без перестановки строк метод мог бы сломаться:
// в первом элементе стоит 0, поэтому нужен выбор ведущего
// элемента и обмен строк.
 
void testPivotingIsReallyUsed()
{
    std::vector<std::vector<double>> A = {
        {0.0, 2.0},
        {1.0, 3.0}
    };
    std::vector<double> b = {4.0, 5.0};

    GaussianSolver solver;
    const GaussResult result = solver.solve(A, b);

    require(result.type == SolutionType::Unique, "Ожидалось единственное решение.");
    require(!result.isDegenerate, "Матрица не должна быть вырожденной.");

    const std::vector<double> expected = {-1.0, 2.0};
    require(vectorApproxEqual(result.x, expected, 1e-9),
            "Неверное решение при необходимости выбора ведущего элемента.");
}

 
// Тест на несовместную систему.
// Пример:
// x + y = 2
// 2x + 2y = 5
//
// Вторая строка противоречит первой, решений нет.
 
void testNoSolution()
{
    std::vector<std::vector<double>> A = {
        {1.0, 1.0},
        {2.0, 2.0}
    };
    std::vector<double> b = {2.0, 5.0};

    GaussianSolver solver;
    const GaussResult result = solver.solve(A, b);

    require(result.type == SolutionType::None, "Ожидалось отсутствие решений.");
    require(result.isDegenerate, "Матрица должна быть вырожденной.");
    require(result.rankA < result.rankAugmented, "Должно быть rank(A) < rank([A|b]).");
}

 
// Тест на бесконечно много решений.
// Пример:
// x + y = 2
// 2x + 2y = 4
//
// Вторая строка линейно зависит от первой.
 
void testInfiniteSolutions()
{
    std::vector<std::vector<double>> A = {
        {1.0, 1.0},
        {2.0, 2.0}
    };
    std::vector<double> b = {2.0, 4.0};

    GaussianSolver solver;
    const GaussResult result = solver.solve(A, b);

    require(result.type == SolutionType::Infinite, "Ожидалось бесконечно много решений.");
    require(result.isDegenerate, "Матрица должна быть вырожденной.");
    require(result.rankA == result.rankAugmented, "Для совместной вырожденной системы ранги должны совпадать.");

    // Проверяем, что найденное частное решение действительно удовлетворяет системе.
    const double residual = residualInfNorm(A, result.x, b);
    require(residual < 1e-9, "Частное решение не удовлетворяет системе.");
}

 
// Тест на единичной матрице.
// Ожидаем, что решение равно правой части.
 
void testIdentityMatrix()
{
    const int n = 5;
    std::vector<std::vector<double>> A = makeIdentityMatrix(n);
    std::vector<double> b = {10.0, -2.0, 3.5, 7.0, 0.25};

    GaussianSolver solver;
    const GaussResult result = solver.solve(A, b);

    require(result.type == SolutionType::Unique, "Ожидалось единственное решение.");
    require(!result.isDegenerate, "Единичная матрица невырождена.");
    require(vectorApproxEqual(result.x, b, 1e-12), "Для единичной матрицы x должно совпадать с b.");
}

 
// Тест на матрице Гильберта.
//
// Важно:
// матрица Гильберта очень плохо обусловлена.
// Поэтому проверять точное совпадение координат решения
// не всегда разумно. Гораздо правильнее проверять невязку.
//
// Здесь мы задаём x_true = [1, 1, ..., 1],
// считаем b = H * x_true,
// затем решаем систему и проверяем:
// 1) что решение найдено,
// 2) что невязка маленькая.
 
void testHilbertMatrix(int n, double maxResidualAllowed)
{
    std::vector<std::vector<double>> H = makeHilbertMatrix(n);
    std::vector<double> xTrue(n, 1.0);
    std::vector<double> b = multiply(H, xTrue);

    GaussianSolver solver(1e-14);
    const GaussResult result = solver.solve(H, b);

    require(result.type == SolutionType::Unique,
            "Для матрицы Гильберта ожидалось единственное решение.");

    const double residual = residualInfNorm(H, result.x, b);
    require(residual < maxResidualAllowed,
            "Слишком большая невязка на матрице Гильберта.");
}

 
// Большой тест на диагональной матрице.
// Такой случай удобен тем, что:
// - матрица точно невырождена,
// - решение известно заранее,
// - тест не такой "ядовитый", как случайная плотная матрица.
 
void testLargeDiagonalMatrix(int n)
{
    std::vector<std::vector<double>> A = makeDiagonalMatrix(n);

    // Возьмём истинное решение x = [1, 1, ..., 1].
    std::vector<double> xTrue(n, 1.0);

    // Тогда b просто равен диагонали.
    std::vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        b[i] = static_cast<double>(i + 1);
    }

    GaussianSolver solver;

    const auto start = std::chrono::steady_clock::now();
    const GaussResult result = solver.solve(A, b);
    const auto finish = std::chrono::steady_clock::now();

    const std::chrono::duration<double> elapsed = finish - start;

    require(result.type == SolutionType::Unique, "Ожидалось единственное решение для большой диагональной матрицы.");
    require(!result.isDegenerate, "Диагональная матрица с ненулевой диагональю невырождена.");

    const double residual = residualInfNorm(A, result.x, b);
    require(residual < 1e-9, "Слишком большая невязка на большой диагональной матрице.");

    std::cout << "       Размер n = " << n
              << ", время = " << elapsed.count() << " сек.\n";
}

 
// Главная функция тестов.
 
int main()
{
    std::cout << "Запуск тестов метода Гаусса\n";
    std::cout << "----------------------------------------\n";

    runTest("Базовый случай 2x2", testBasicUnique);
    runTest("Проверка выбора ведущего элемента", testPivotingIsReallyUsed);
    runTest("Несовместная система", testNoSolution);
    runTest("Бесконечно много решений", testInfiniteSolutions);
    runTest("Единичная матрица", testIdentityMatrix);

    // Матрицы Гильберта.
    runTest("Матрица Гильберта 5x5", []() {
        testHilbertMatrix(5, 1e-10);
    });

    runTest("Матрица Гильберта 8x8", []() {
        testHilbertMatrix(8, 1e-8);
    });

    runTest("Матрица Гильберта 13x13", []() {
        testHilbertMatrix(13, 1e-3);
    });

    runTest("Матрица Гильберта 20x20", []() {
        testHilbertMatrix(20, 1e-3);
    });

    runTest("Матрица Гильберта 30x30", []() {
        testHilbertMatrix(30, 1e-3);
    });

    // Большой тест 1000x1000.
    // Он уже вполне серьёзный по памяти и времени,
    // но обычно ещё выполним на обычной машине.
    runTest("Большая диагональная матрица 1000x1000", []() {
        testLargeDiagonalMatrix(1000);
    });

    // Очень большой тест 10000x10000.
    //
    // ВАЖНО:
    // Для плотного метода Гаусса такой размер уже крайне тяжёлый:
    // - память очень большая,
    // - время у классического метода O(n^3).
    //
    // Поэтому этот тест по умолчанию отключён.
    // Его можно включить при компиляции макросом:
    // -DRUN_VERY_LARGE_TESTS
    //
    // Это честнее, чем делать вид, что такой тест "лёгкий".
    // На практике для 10000x10000 плотный Гаусс - это уже
    // очень тяжёлая задача.
#ifdef RUN_VERY_LARGE_TESTS
    runTest("Очень большая диагональная матрица 10000x10000", []() {
        testLargeDiagonalMatrix(10000);
    });
#else
    std::cout << "[SKIP] Очень большая диагональная матрица 10000x10000 "
                 "(отключено по умолчанию, включается через -DRUN_VERY_LARGE_TESTS)\n";
#endif

    std::cout << "----------------------------------------\n";
    std::cout << "Пройдено тестов: " << g_passedTests
              << " из " << g_totalTests << "\n";

    if (g_passedTests == g_totalTests) {
        std::cout << "Все тесты успешно пройдены.\n";
        return 0;
    }

    std::cout << "Есть непройденные тесты.\n";
    return 1;
}
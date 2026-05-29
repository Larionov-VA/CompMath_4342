#include "jacobi_solver.hpp"

#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using namespace jacobi;

namespace {

void printHeader(const std::string& title) {
    std::cout << "\n==============================\n";
    std::cout << title << "\n";
    std::cout << "==============================\n";
}

void assertClose(double actual, double expected, double eps, const std::string& message) {
    if (std::fabs(actual - expected) > eps) {
        std::cerr << "Провал проверки: " << message
                  << ", actual = " << actual
                  << ", expected = " << expected << "\n";
        std::abort();
    }
}

void assertVectorClose(const Vector& actual, const Vector& expected, double eps, const std::string& message) {
    if (actual.size() != expected.size()) {
        std::cerr << "Провал проверки: размеры векторов различаются. " << message << "\n";
        std::abort();
    }
    for (std::size_t i = 0; i < actual.size(); ++i) {
        if (std::fabs(actual[i] - expected[i]) > eps) {
            std::cerr << "Провал проверки: " << message
                      << ", index = " << i
                      << ", actual = " << actual[i]
                      << ", expected = " << expected[i] << "\n";
            std::abort();
        }
    }
}

void testBasicCase() {
    printHeader("Тест 1. Базовый диагонально доминирующий случай");

    Matrix a = {
        {10.0, -1.0,  2.0},
        {-1.0, 11.0, -1.0},
        { 2.0, -1.0, 10.0}
    };
    Vector exact = {1.0, 2.0, -1.0};
    Vector b = {6.0, 22.0, -10.0};

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(a, b, 1e-12, 10000, true);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assert(result.guaranteed_convergence);
    assertVectorClose(result.x, exact, 1e-8, "Базовый тест: неверное решение");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

void testNoSolution() {
    printHeader("Тест 2. Несовместная система");

    Matrix a = {
        {1.0, 1.0},
        {2.0, 2.0}
    };
    Vector b = {2.0, 5.0};

    const SystemStatus status = analyzeSystem(a, b);
    assert(status.type == SolutionType::NONE);

    const JacobiResult result = solveJacobi(a, b);
    assert(!result.converged);

    std::cout << "OK, программа правильно определила отсутствие решения.\n";
}

void testInfiniteSolutions() {
    printHeader("Тест 3. Система с бесконечно многими решениями");

    Matrix a = {
        {1.0, 1.0},
        {2.0, 2.0}
    };
    Vector b = {2.0, 4.0};

    const SystemStatus status = analyzeSystem(a, b);
    assert(status.type == SolutionType::INFINITE);

    const JacobiResult result = solveJacobi(a, b);
    assert(!result.converged);

    std::cout << "OK, программа правильно определила бесконечное число решений.\n";
}

void testHilbert2() {
    printHeader("Тест 4. Матрица Гильберта 2x2");

    Matrix a;
    Vector b;
    Vector exact = {1.0, 1.0};
    generateHilbertSystem(2, a, b, exact);

    const JacobiResult result = solveJacobi(a, b, 1e-12, 100000, true);

    assert(result.converged);
    assertVectorClose(result.x, exact, 1e-6, "Матрица Гильберта 2x2: неверное решение");
    std::cout << "OK, для H2 метод Якоби сошёлся.\n";
}

void testHilbert3() {
    printHeader("Тест 5. Матрица Гильберта 3x3");

    Matrix a;
    Vector b;
    Vector exact = {1.0, 1.0, 1.0};
    generateHilbertSystem(3, a, b, exact);

    const JacobiResult result = solveJacobi(a, b, 1e-12, 500, true);

    assert(!result.guaranteed_convergence);
    std::cout << "Сообщение программы: " << result.message << "\n";
    std::cout << "OK, программа корректно показала, что гарантии сходимости нет.\n";
}

void testDenseDD() {
    printHeader("Тест 6. Плотная диагонально доминирующая матрица 150x150");

    const int n = 150;
    Matrix a;
    Vector b;
    Vector exact(n, 1.0);
    generateDenseDiagonallyDominantSystem(n, a, b, exact, 42u);

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(a, b, 1e-10, 10000, true);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assert(result.guaranteed_convergence);
    assertClose(result.residual_norm, 0.0, 1e-6, "Плотная диагонально доминирующая система: большая невязка");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

void testLargeTridiag1000() {
    printHeader("Тест 7. Большая система 1000x1000 (трёхдиагональная)");

    const int n = 1000;
    Vector exact(n, 1.0);
    TridiagonalSystem s = generateLargeTridiagonalSystem(n, exact);

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(s, 1e-10, 100000);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assertClose(result.residual_norm, 0.0, 1e-6, "Большая трёхдиагональная система: большая невязка");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

void testLargePentadiag2000() {
    printHeader("Тест 8. Большая система 2000x2000 (пятидиагональная)");

    const int n = 2000;
    Vector exact(n, 1.0);
    BandedSystem s = generatePentadiagonalSystem(n, exact, 7u);

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(s, 1e-10, 100000);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assert(result.guaranteed_convergence);
    assertClose(result.residual_norm, 0.0, 1e-6, "Пятидиагональная система: большая невязка");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

void testLargeBanded3000() {
    printHeader("Тест 9. Большая система 3000x3000 (ленточная, ширина 3 снизу и 4 сверху)");

    const int n = 3000;
    Vector exact(n, 1.0);
    BandedSystem s = generateBandedDiagonallyDominantSystem(n, 3, 4, exact, 11u);

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(s, 1e-10, 100000);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assert(result.guaranteed_convergence);
    assertClose(result.residual_norm, 0.0, 1e-6, "Ленточная система: большая невязка");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

void testLargeSymmBanded10000() {
    printHeader("Тест 10. Большая система 10000x10000 (симметричная SPD-ленточная)");

    const int n = 10000;
    Vector exact(n, 1.0);
    BandedSystem s = generateSymmetricPositiveDefiniteBandedSystem(n, 3, exact, 13u);

    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = solveJacobi(s, 1e-8, 100000);
    const auto finish = std::chrono::steady_clock::now();

    assert(result.converged);
    assert(result.guaranteed_convergence);
    assertClose(result.residual_norm, 0.0, 1e-5, "Симметричная SPD-ленточная система: большая невязка");

    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    std::cout << "OK, итераций = " << result.iterations
              << ", невязка = " << result.residual_norm
              << ", время = " << ms << " мс\n";
}

} // namespace

int main() {
    std::cout << std::fixed << std::setprecision(10);

    testBasicCase();
    testNoSolution();
    testInfiniteSolutions();
    testHilbert2();
    testHilbert3();
    testDenseDD();
    testLargeTridiag1000();
    testLargePentadiag2000();
    testLargeBanded3000();
    testLargeSymmBanded10000();

    std::cout << "\nВсе тесты успешно пройдены.\n";
    return 0;
}

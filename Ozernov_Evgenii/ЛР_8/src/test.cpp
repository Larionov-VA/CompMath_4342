#include "jacobi_solver.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Удобная структура для описания одного теста.
struct TestCase {
    int id = 0;
    std::string name;
    SparseMatrix a;
    std::vector<double> b;
    std::vector<double> referenceX;
    JacobiConfig config;
    std::string inputDescription;

    // Проверка результата теста.
    // Возвращает true/false и пишет пояснение в outReason.
    std::function<bool(const SolveResult&, std::string& outReason)> check;
};

// Умножение разреженной матрицы на вектор (для генерации эталонного b).
std::vector<double> multiply(const SparseMatrix& a, const std::vector<double>& x) {
    const std::size_t n = a.size();
    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (const auto& [col, value] : a.row(i)) {
            sum += value * x[col];
        }
        result[i] = sum;
    }
    return result;
}

// Максимальная абсолютная ошибка между векторами.
double maxAbsError(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double mx = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        mx = std::max(mx, std::abs(a[i] - b[i]));
    }
    return mx;
}

// Форматирование первых элементов вектора для отчёта.
std::string shortVector(const std::vector<double>& v, std::size_t maxItems = 12) {
    std::ostringstream oss;
    oss << "[";
    for (std::size_t i = 0; i < v.size() && i < maxItems; ++i) {
        if (i > 0) {
            oss << ", ";
        }
        oss << std::setprecision(10) << v[i];
    }
    if (v.size() > maxItems) {
        oss << ", ...";
    }
    oss << "]";
    return oss.str();
}

// Генератор циклической ленточной матрицы:
// A[i][i] = diagValue, а также несколько "смещённых" диагоналей по модулю n.
// Такой тип матриц большой, разреженный и не сводится к простой диагонали.
SparseMatrix buildCyclicBandedMatrix(
    std::size_t n,
    double diagValue,
    const std::vector<std::pair<int, double>>& offsetsAndValues) {
    SparseMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) {
        m.set(i, i, diagValue);
        for (const auto& [offset, value] : offsetsAndValues) {
            const long long jSigned =
                (static_cast<long long>(i) + offset + static_cast<long long>(n)) %
                static_cast<long long>(n);
            const std::size_t j = static_cast<std::size_t>(jSigned);
            m.set(i, j, value);
        }
    }
    return m;
}

// Генератор пятидиагональной матрицы (две соседние и две "вторые" диагонали).
SparseMatrix buildPentadiagonalMatrix(
    std::size_t n,
    double diagValue,
    double firstOffDiag,
    double secondOffDiag) {
    SparseMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) {
        m.set(i, i, diagValue);
        if (i >= 1) {
            m.set(i, i - 1, firstOffDiag);
        }
        if (i + 1 < n) {
            m.set(i, i + 1, firstOffDiag);
        }
        if (i >= 2) {
            m.set(i, i - 2, secondOffDiag);
        }
        if (i + 2 < n) {
            m.set(i, i + 2, secondOffDiag);
        }
    }
    return m;
}

// Генератор разреженной матрицы с "дальними" связями:
// для каждой строки добавляются связи с несколькими удалёнными столбцами,
// что имитирует более сложную структуру, чем локальные диагонали.
SparseMatrix buildFarCoupledSparseMatrix(std::size_t n, double diagValue) {
    SparseMatrix m(n);
    for (std::size_t i = 0; i < n; ++i) {
        // Набор детерминированных удалённых смещений.
        // Значения подобраны так, чтобы сохранялась диагональная доминантность.
        const std::size_t j1 = (i + 17) % n;
        const std::size_t j2 = (i + 103) % n;
        const std::size_t j3 = (i + 509) % n;
        const std::size_t j4 = (i + n - 31) % n;
        const std::size_t j5 = (i + n - 211) % n;

        m.set(i, j1, -0.8);
        m.set(i, j2, 0.6);
        m.set(i, j3, -0.4);
        m.set(i, j4, 0.7);
        m.set(i, j5, -0.5);

        // Диагональ делаем достаточно большой для устойчивой сходимости Якоби.
        m.set(i, i, diagValue);
    }
    return m;
}

// Создание базового конфигурационного шаблона метода Якоби.
JacobiConfig baseConfig() {
    JacobiConfig cfg;
    cfg.tolerance = 1e-10;
    cfg.maxIterations = 200000;
    cfg.zeroEps = 1e-12;
    cfg.exactRankCheckLimit = 250;
    cfg.pivotStrategy = PivotStrategy::PARTIAL_COLUMN_PIVOTING;
    return cfg;
}

} // namespace

int main() {
    std::vector<TestCase> tests;
    tests.reserve(10);

    // Тест 1: базовый случай с единственным решением.
    {
        SparseMatrix a(3);
        a.set(0, 0, 10.0); a.set(0, 1, 1.0);  a.set(0, 2, 1.0);
        a.set(1, 0, 2.0);  a.set(1, 1, 10.0); a.set(1, 2, 1.0);
        a.set(2, 0, 2.0);  a.set(2, 1, 2.0);  a.set(2, 2, 10.0);

        std::vector<double> xTrue = {1.0, 1.0, 1.0};
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.tolerance = 1e-12;

        tests.push_back(TestCase{
            1,
            "Базовый случай: единственное решение (3x3)",
            a,
            b,
            xTrue,
            cfg,
            "n=3, плотная диагонально-доминантная матрица, x_true=[1,1,1]",
            [](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Ожидалась успешная сходимость, но её нет.";
                    return false;
                }
                const double err = maxAbsError(res.x, std::vector<double>{1.0, 1.0, 1.0});
                if (err > 1e-8) {
                    reason = "Слишком большая ошибка решения: " + std::to_string(err);
                    return false;
                }
                reason = "Решение найдено, ошибка в пределах нормы.";
                return true;
            }
        });
    }

    // Тест 2: система без решений (несовместная).
    {
        SparseMatrix a(2);
        a.set(0, 0, 1.0); a.set(0, 1, 1.0);
        a.set(1, 0, 2.0); a.set(1, 1, 2.0);
        std::vector<double> b = {2.0, 5.0};
        JacobiConfig cfg = baseConfig();

        tests.push_back(TestCase{
            2,
            "Базовый случай: решений нет",
            a,
            b,
            {},
            cfg,
            "n=2, зависимые строки A, но несовместный b",
            [](const SolveResult& res, std::string& reason) {
                if (res.solutionType != SolutionType::NONE) {
                    reason = "Ожидался тип решения NONE (нет решений).";
                    return false;
                }
                reason = "Несовместность корректно обнаружена.";
                return true;
            }
        });
    }

    // Тест 3: бесконечно много решений.
    {
        SparseMatrix a(2);
        a.set(0, 0, 1.0); a.set(0, 1, 1.0);
        a.set(1, 0, 2.0); a.set(1, 1, 2.0);
        std::vector<double> b = {2.0, 4.0};
        JacobiConfig cfg = baseConfig();

        tests.push_back(TestCase{
            3,
            "Базовый случай: бесконечно много решений",
            a,
            b,
            {},
            cfg,
            "n=2, зависимые строки A и согласованный b",
            [](const SolveResult& res, std::string& reason) {
                if (res.solutionType != SolutionType::INFINITE) {
                    reason = "Ожидался тип решения INFINITE.";
                    return false;
                }
                reason = "Случай бесконечного множества решений корректно обнаружен.";
                return true;
            }
        });
    }

    // Тест 4: "абсурдно большая" матрица 1000x1000.
    // Используется трёхдиагональная матрица, чтобы тест был реалистичен по памяти.
    {
        const std::size_t n = 1000;
        SparseMatrix a = buildTridiagonalMatrix(n, 10.0, -1.0);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;  // Для этой матрицы перестановки не нужны.
        cfg.tolerance = 1e-8;
        cfg.maxIterations = 5000;

        tests.push_back(TestCase{
            4,
            "Большая матрица: 1000x1000 (трёхдиагональная)",
            a,
            b,
            xTrue,
            cfg,
            "n=1000, tridiagonal(diag=10, offdiag=-1), x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Для 1000x1000 ожидалась сходимость, но она не достигнута.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 1e-5) {
                    reason = "Слишком большая ошибка на 1000x1000: " + std::to_string(err);
                    return false;
                }
                reason = "Сходимость и точность для 1000x1000 подтверждены.";
                return true;
            }
        });
    }

    // Тест 5: "абсурдно большая" матрица 10000x10000.
    {
        const std::size_t n = 10000;
        SparseMatrix a = buildTridiagonalMatrix(n, 10.0, -1.0);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;
        cfg.tolerance = 1e-6;
        cfg.maxIterations = 2000;
        cfg.exactRankCheckLimit = 150; // Явно отключаем дорогой ранговый анализ для такого n.

        tests.push_back(TestCase{
            5,
            "Большая матрица: 10000x10000 (трёхдиагональная)",
            a,
            b,
            xTrue,
            cfg,
            "n=10000, tridiagonal(diag=10, offdiag=-1), x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Для 10000x10000 ожидалась сходимость, но её нет.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 1e-4) {
                    reason = "Слишком большая ошибка на 10000x10000: " + std::to_string(err);
                    return false;
                }
                reason = "Сходимость и точность для 10000x10000 подтверждены.";
                return true;
            }
        });
    }

    // Тест 6: матрица Гильберта (обычно сложна для Якоби, часто даёт расходимость).
    {
        const std::size_t n = 8;
        SparseMatrix a = buildHilbertMatrix(n);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.tolerance = 1e-12;
        cfg.maxIterations = 3000;

        tests.push_back(TestCase{
            6,
            "Специальная матрица: Гильберта 8x8",
            a,
            b,
            xTrue,
            cfg,
            "n=8, Hilbert matrix, x_true=ones",
            [](const SolveResult& res, std::string& reason) {
                // Для классического Якоби матрица Гильберта, как правило, неблагоприятна.
                // Здесь считаем корректным исходом отсутствие сходимости.
                if (res.converged) {
                    reason = "Неожиданно получена сходимость на Hilbert 8x8.";
                    return false;
                }
                reason = "Ожидаемая трудность/расходимость на матрице Гильберта подтверждена.";
                return true;
            }
        });
    }

    // Тест 7: большая циклическая ленточная матрица 4000x4000.
    {
        const std::size_t n = 4000;
        SparseMatrix a = buildCyclicBandedMatrix(
            n,
            10.0,
            {
                {1, -1.0}, {-1, -1.0},
                {17, 0.4}, {-17, 0.4},
                {101, -0.3}, {-101, -0.3}
            }
        );
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;
        cfg.tolerance = 1e-8;
        cfg.maxIterations = 4000;
        cfg.exactRankCheckLimit = 120;

        tests.push_back(TestCase{
            7,
            "Большая матрица: 4000x4000 (циклическая ленточная, 7 диагоналей)",
            a,
            b,
            xTrue,
            cfg,
            "n=4000, cyclic banded offsets={+-1, +-17, +-101}, diag=10, x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Ожидалась сходимость на циклической ленточной матрице 4000x4000.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 2e-5) {
                    reason = "Слишком большая ошибка на cyclic-banded 4000x4000: " + std::to_string(err);
                    return false;
                }
                reason = "Циклическая ленточная матрица 4000x4000 успешно пройдена.";
                return true;
            }
        });
    }

    // Тест 8: очень большая циклическая ленточная матрица 12000x12000.
    {
        const std::size_t n = 12000;
        SparseMatrix a = buildCyclicBandedMatrix(
            n,
            11.0,
            {
                {1, -1.0}, {-1, -1.0},
                {37, 0.5}, {-37, 0.5}
            }
        );
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;
        cfg.tolerance = 2e-6;
        cfg.maxIterations = 2500;
        cfg.exactRankCheckLimit = 80;

        tests.push_back(TestCase{
            8,
            "Большая матрица: 12000x12000 (циклическая ленточная, 5 диагоналей)",
            a,
            b,
            xTrue,
            cfg,
            "n=12000, cyclic banded offsets={+-1, +-37}, diag=11, x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Ожидалась сходимость на циклической ленточной матрице 12000x12000.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 3e-4) {
                    reason = "Слишком большая ошибка на cyclic-banded 12000x12000: " + std::to_string(err);
                    return false;
                }
                reason = "Циклическая ленточная матрица 12000x12000 успешно пройдена.";
                return true;
            }
        });
    }

    // Тест 9: большая пятидиагональная матрица 6000x6000.
    {
        const std::size_t n = 6000;
        SparseMatrix a = buildPentadiagonalMatrix(n, 12.0, -1.0, 0.7);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;
        cfg.tolerance = 1e-7;
        cfg.maxIterations = 3500;
        cfg.exactRankCheckLimit = 100;

        tests.push_back(TestCase{
            9,
            "Большая матрица: 6000x6000 (пятидиагональная)",
            a,
            b,
            xTrue,
            cfg,
            "n=6000, pentadiagonal(diag=12, off1=-1, off2=0.7), x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Ожидалась сходимость на пятидиагональной матрице 6000x6000.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 1e-4) {
                    reason = "Слишком большая ошибка на pentadiagonal 6000x6000: " + std::to_string(err);
                    return false;
                }
                reason = "Пятидиагональная матрица 6000x6000 успешно пройдена.";
                return true;
            }
        });
    }

    // Тест 10: большая разреженная матрица с дальними связями 9000x9000.
    {
        const std::size_t n = 9000;
        SparseMatrix a = buildFarCoupledSparseMatrix(n, 14.0);
        std::vector<double> xTrue(n, 1.0);
        std::vector<double> b = multiply(a, xTrue);
        JacobiConfig cfg = baseConfig();
        cfg.pivotStrategy = PivotStrategy::NONE;
        cfg.tolerance = 2e-6;
        cfg.maxIterations = 2500;
        cfg.exactRankCheckLimit = 80;

        tests.push_back(TestCase{
            10,
            "Большая матрица: 9000x9000 (разреженная с дальними связями)",
            a,
            b,
            xTrue,
            cfg,
            "n=9000, sparse far-coupled offsets={+17,+103,+509,-31,-211}, diag=14, x_true=ones",
            [xTrue](const SolveResult& res, std::string& reason) {
                if (!res.success || !res.converged) {
                    reason = "Ожидалась сходимость на far-coupled матрице 9000x9000.";
                    return false;
                }
                const double err = maxAbsError(res.x, xTrue);
                if (err > 4e-4) {
                    reason = "Слишком большая ошибка на far-coupled 9000x9000: " + std::to_string(err);
                    return false;
                }
                reason = "Разреженная матрица с дальними связями 9000x9000 успешно пройдена.";
                return true;
            }
        });
    }

    std::ofstream report("test_report.txt", std::ios::trunc);
    if (!report.is_open()) {
        std::cerr << "Не удалось открыть test_report.txt для записи.\n";
        return 1;
    }

    report << "Отчёт по тестированию метода Якоби\n";
    report << "===================================\n\n";

    int passed = 0;
    for (const TestCase& test : tests) {
        const auto start = std::chrono::high_resolution_clock::now();
        const SolveResult res = solveByJacobi(test.a, test.b, test.config);
        const auto finish = std::chrono::high_resolution_clock::now();
        const auto elapsedMs =
            std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

        std::string reason;
        const bool ok = test.check(res, reason);
        if (ok) {
            ++passed;
        }

        std::cout << "Test " << test.id << " - " << test.name << ": " << (ok ? "OK" : "Fail") << "\n";

        report << "Тест #" << test.id << ": " << test.name << "\n";
        report << "Состояние: " << (ok ? "OK" : "Fail") << "\n";
        report << "Комментарий проверки: " << reason << "\n\n";

        report << "Ввод:\n";
        report << test.inputDescription << "\n";
        report << "Параметры: tolerance=" << test.config.tolerance
               << ", maxIterations=" << test.config.maxIterations
               << ", pivotStrategy="
               << (test.config.pivotStrategy == PivotStrategy::PARTIAL_COLUMN_PIVOTING ? "PARTIAL" : "NONE")
               << ", zeroEps=" << test.config.zeroEps
               << ", exactRankCheckLimit=" << test.config.exactRankCheckLimit << "\n\n";

        report << "Вывод:\n";
        report << "success=" << (res.success ? "true" : "false")
               << ", converged=" << (res.converged ? "true" : "false")
               << ", iterations=" << res.iterations
               << ", achievedDelta=" << std::setprecision(16) << res.achievedDelta << "\n";
        report << "estimatedMemoryBytes=" << res.estimatedMemoryBytes << "\n";
        report << "message=" << res.message << "\n";
        if (!res.x.empty()) {
            report << "x(first values)=" << shortVector(res.x, 10) << "\n";
        }
        report << "\n";

        report << "Время (мс): " << elapsedMs << "\n";
        report << "Память (байт): " << res.estimatedMemoryBytes << "\n";
        report << "-----------------------------------\n\n";
    }

    std::cout << "Итог: " << passed << "/" << static_cast<int>(tests.size()) << " тестов успешно.\n";
    std::cout << "Подробный отчёт: test_report.txt\n";

    return (passed == static_cast<int>(tests.size())) ? 0 : 1;
}

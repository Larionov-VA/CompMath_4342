#include "gauss_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#define PSAPI_VERSION 2
#include <windows.h>
#include <psapi.h>
#endif

namespace {

struct TestRunResult {
    bool ok = true;
    std::string input;
    std::string output;
    std::string error;
    double elapsed_ms = 0.0;             // Время выполнения теста.
    std::uint64_t working_set_bytes = 0; // Рабочий набор после теста.
    std::uint64_t peak_delta_bytes = 0;  // Прирост пикового рабочего набора в тесте.
};

struct MemorySnapshot {
    std::uint64_t working_set_bytes = 0;
    std::uint64_t peak_working_set_bytes = 0;
};

// Сравнение чисел с учетом допуска eps.
bool approxEqual(double a, double b, double eps = 1e-7) {
    return std::fabs(a - b) <= eps;
}

// Проверка условия внутри теста: если условие ложно,
// тест помечается как неуспешный и сохраняется текст ошибки.
void expect(bool condition, const std::string& message, TestRunResult& run) {
    if (!condition && run.ok) {
        run.ok = false;
        run.error = message;
    }
}

// Преобразование типа результата решателя в читаемый русский текст.
std::string solutionTypeToText(SolutionType type) {
    switch (type) {
        case SolutionType::Unique:
            return "Единственное решение";
        case SolutionType::Infinite:
            return "Бесконечно много решений";
        case SolutionType::NoSolution:
            return "Решений нет";
        case SolutionType::InvalidInput:
            return "Некорректный ввод";
    }
    return "Неизвестный статус";
}

// Чтение текущего и пикового потребления памяти процесса.
MemorySnapshot getMemorySnapshot() {
    MemorySnapshot snapshot;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters{};
    if (K32GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)) != 0) {
        snapshot.working_set_bytes = static_cast<std::uint64_t>(counters.WorkingSetSize);
        snapshot.peak_working_set_bytes = static_cast<std::uint64_t>(counters.PeakWorkingSetSize);
    }
#endif
    return snapshot;
}

// Формат размера памяти в удобном виде.
std::string formatBytes(std::uint64_t bytes) {
    static const char* units[] = {"Б", "КБ", "МБ", "ГБ"};
    long double value = static_cast<long double>(bytes);
    int unit = 0;

    while (value >= 1024.0L && unit < 3) {
        value /= 1024.0L;
        ++unit;
    }

    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << static_cast<double>(value) << " " << units[unit];
    return out.str();
}

// Краткая печать вектора: для больших размеров показываем только начало и конец.
std::string formatVectorBrief(const std::vector<double>& v, std::size_t edge = 3) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    out << "[";
    if (v.size() <= 2 * edge) {
        for (std::size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                out << ", ";
            }
            out << v[i];
        }
    } else {
        for (std::size_t i = 0; i < edge; ++i) {
            if (i > 0) {
                out << ", ";
            }
            out << v[i];
        }

        out << ", ..., ";

        for (std::size_t i = v.size() - edge; i < v.size(); ++i) {
            if (i > v.size() - edge) {
                out << ", ";
            }
            out << v[i];
        }
    }
    out << "]";

    return out.str();
}

// Краткая печать одной строки матрицы: показываем только начало и конец столбцов.
std::string formatRowBrief(const std::vector<double>& row, std::size_t edge = 3) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    out << "[";
    if (row.size() <= 2 * edge) {
        for (std::size_t j = 0; j < row.size(); ++j) {
            if (j > 0) {
                out << ", ";
            }
            out << row[j];
        }
    } else {
        for (std::size_t j = 0; j < edge; ++j) {
            if (j > 0) {
                out << ", ";
            }
            out << row[j];
        }

        out << ", ..., ";

        for (std::size_t j = row.size() - edge; j < row.size(); ++j) {
            if (j > row.size() - edge) {
                out << ", ";
            }
            out << row[j];
        }
    }
    out << "]";

    return out.str();
}

// Краткая печать матрицы: для больших размеров выводим только первые и последние строки.
std::string formatMatrixBrief(const std::vector<std::vector<double>>& A, std::size_t row_edge = 2, std::size_t col_edge = 3) {
    std::ostringstream out;
    const std::size_t m = A.size();
    const std::size_t n = A.empty() ? 0 : A[0].size();

    out << "Размер матрицы A: " << m << "x" << n << "\n";

    if (m == 0 || n == 0) {
        out << "[]\n";
        return out.str();
    }

    if (m <= 2 * row_edge) {
        for (std::size_t i = 0; i < m; ++i) {
            out << "  строка " << i << ": " << formatRowBrief(A[i], col_edge) << "\n";
        }
    } else {
        for (std::size_t i = 0; i < row_edge; ++i) {
            out << "  строка " << i << ": " << formatRowBrief(A[i], col_edge) << "\n";
        }

        out << "  ...\n";

        for (std::size_t i = m - row_edge; i < m; ++i) {
            out << "  строка " << i << ": " << formatRowBrief(A[i], col_edge) << "\n";
        }
    }

    return out.str();
}

// Формирование текстового блока с входными данными теста.
std::string formatSystemBrief(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    std::ostringstream out;
    out << formatMatrixBrief(A);
    out << "Вектор b: " << formatVectorBrief(b) << "\n";
    return out.str();
}

// Формирование текстового блока с ответом решателя.
std::string formatSolverOutput(const GaussResult& result) {
    std::ostringstream out;

    out << "Статус: " << solutionTypeToText(result.type) << "\n";
    out << "Сообщение: " << result.message << "\n";
    out << "Ранг: " << result.rank << "\n";
    out << "Невырожденность: " << (result.nonsingular ? "да" : "нет") << "\n";

    if (result.type == SolutionType::Unique || result.type == SolutionType::Infinite) {
        out << "Вектор x: " << formatVectorBrief(result.x) << "\n";
    }

    return out.str();
}

// Генератор матрицы Гильберта размерности n.
std::vector<std::vector<double>> makeHilbert(std::size_t n) {
    std::vector<std::vector<double>> h(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            h[i][j] = 1.0 / static_cast<double>(i + j + 1);
        }
    }
    return h;
}

// Умножение матрицы на вектор: b = A * x.
std::vector<double> matVec(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    const std::size_t m = A.size();
    const std::size_t n = x.size();
    std::vector<double> b(m, 0.0);

    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }

    return b;
}

// Базовый тест: квадратная невырожденная система с единственным решением.
TestRunResult testUniqueBasic() {
    TestRunResult run;

    std::vector<std::vector<double>> A{{2.0, 1.0}, {5.0, 7.0}};
    std::vector<double> b{11.0, 13.0};

    run.input = formatSystemBrief(A, b);

    const auto result = solveGaussian(A, b);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::Unique, "Ожидалось единственное решение", run);
    expect(result.nonsingular, "Матрица должна быть невырожденной", run);
    expect(approxEqual(result.x[0], 64.0 / 9.0, 1e-8), "Неверное значение x0", run);
    expect(approxEqual(result.x[1], -29.0 / 9.0, 1e-8), "Неверное значение x1", run);

    return run;
}

// Тест на несовместную систему (решений нет).
TestRunResult testNoSolution() {
    TestRunResult run;

    std::vector<std::vector<double>> A{{1.0, 1.0}, {1.0, 1.0}};
    std::vector<double> b{1.0, 2.0};

    run.input = formatSystemBrief(A, b);

    const auto result = solveGaussian(A, b);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::NoSolution, "Ожидался статус «решений нет»", run);
    expect(!result.nonsingular, "Невырожденность невозможна при несовместной системе", run);

    return run;
}

// Тест на бесконечно много решений.
TestRunResult testInfiniteSolutions() {
    TestRunResult run;

    std::vector<std::vector<double>> A{{1.0, 1.0}, {2.0, 2.0}};
    std::vector<double> b{2.0, 4.0};

    run.input = formatSystemBrief(A, b);

    const auto result = solveGaussian(A, b);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::Infinite, "Ожидался статус «бесконечно много решений»", run);
    expect(!result.nonsingular, "Матрица должна быть вырожденной", run);

    return run;
}

// Тест на специальной матрице Гильберта (численно сложный пример).
TestRunResult testHilbert() {
    TestRunResult run;

    const std::size_t n = 8;
    auto H = makeHilbert(n);
    std::vector<double> x_true(n, 1.0);
    auto b = matVec(H, x_true);

    run.input = formatSystemBrief(H, b);

    const auto result = solveGaussian(H, b);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::Unique, "Ожидался статус «единственное решение»", run);

    double max_abs_err = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        max_abs_err = std::max(max_abs_err, std::fabs(result.x[i] - x_true[i]));
    }
    expect(max_abs_err < 1e-4, "Слишком большая ошибка восстановления решения", run);

    return run;
}

// Большой тест 1000x1000.
TestRunResult testLarge1000() {
    TestRunResult run;

    const std::size_t n = 1000;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> b(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double d = static_cast<double>(i + 1);
        A[i][i] = d;
        b[i] = d; // Истинное решение: x_i = 1.
    }

    run.input = formatSystemBrief(A, b);

    const auto result = solveGaussian(std::move(A), std::move(b), PivotStrategy::PartialColumnMax, 1e-12, 10000);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::Unique, "Ожидался статус «единственное решение»", run);
    expect(result.nonsingular, "Диагональная матрица с ненулевой диагональю должна быть невырожденной", run);

    for (std::size_t i = 0; i < n; ++i) {
        expect(approxEqual(result.x[i], 1.0, 1e-8), "Неверная компонента решения", run);
    }

    return run;
}

// Большой тест 10000x10000.
TestRunResult testLarge10000() {
    TestRunResult run;

    const std::size_t n = 10000;
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));
    std::vector<double> b(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        A[i][i] = 1.0;
        b[i] = 1.0; // Истинное решение: x_i = 1.
    }

    run.input = formatSystemBrief(A, b);

    const auto result = solveGaussian(std::move(A), std::move(b), PivotStrategy::PartialColumnMax, 1e-12, 10000);
    run.output = formatSolverOutput(result);

    expect(result.type == SolutionType::Unique, "Ожидался статус «единственное решение»", run);
    expect(result.nonsingular, "Единичная матрица должна быть невырожденной", run);

    for (std::size_t i = 0; i < n; ++i) {
        expect(approxEqual(result.x[i], 1.0, 1e-8), "Неверная компонента решения", run);
    }

    return run;
}

// Запись информации по одному тесту в отчетный файл.
void writeTestReport(
    std::ofstream& report,
    std::size_t index,
    const std::string& name,
    const TestRunResult& run
) {
    report << "============================================================\n";
    report << "Тест №" << index << ": " << name << "\n";
    report << "Входные данные:\n";
    report << run.input;
    report << "Ответ решателя:\n";
    report << run.output;
    report << "Время: " << std::fixed << std::setprecision(3) << run.elapsed_ms << " мс\n";
    report << "Память (рабочий набор после теста): " << formatBytes(run.working_set_bytes) << "\n";
    report << "Память (прирост пикового рабочего набора): " << formatBytes(run.peak_delta_bytes) << "\n";
    report << "Итог теста: " << (run.ok ? "OK" : "FAIL") << "\n";
    if (!run.ok) {
        report << "Причина: " << run.error << "\n";
    }
    report << "\n";
}

} // namespace

int main() {
    struct TestCase {
        const char* name;
        TestRunResult (*fn)();
    };

    const std::vector<TestCase> tests{
        {"Базовый уникальный случай", testUniqueBasic},
        {"Несовместная система", testNoSolution},
        {"Бесконечно много решений", testInfiniteSolutions},
        {"Матрица Гильберта", testHilbert},
        {"Большой тест 1000x1000", testLarge1000},
        {"Большой тест 10000x10000", testLarge10000}
    };

    std::ofstream report("test_report.txt", std::ios::out | std::ios::trunc);
    if (!report.is_open()) {
        std::cerr << "Не удалось создать файл отчета test_report.txt\n";
        return 1;
    }

    std::size_t passed = 0;

    for (std::size_t i = 0; i < tests.size(); ++i) {
        const auto& test = tests[i];

        const MemorySnapshot mem_before = getMemorySnapshot();
        const auto time_begin = std::chrono::steady_clock::now();
        TestRunResult run = test.fn();
        const auto time_end = std::chrono::steady_clock::now();
        const MemorySnapshot mem_after = getMemorySnapshot();

        run.elapsed_ms = std::chrono::duration<double, std::milli>(time_end - time_begin).count();
        run.working_set_bytes = mem_after.working_set_bytes;
        run.peak_delta_bytes = (mem_after.peak_working_set_bytes > mem_before.peak_working_set_bytes)
            ? (mem_after.peak_working_set_bytes - mem_before.peak_working_set_bytes)
            : 0;

        writeTestReport(report, i + 1, test.name, run);

        if (run.ok) {
            ++passed;
            std::cout << test.name << ": OK"
                      << " | время: " << std::fixed << std::setprecision(3) << run.elapsed_ms << " мс"
                      << " | память(пик+): " << formatBytes(run.peak_delta_bytes) << "\n";
        } else {
            std::cout << test.name << ": FAIL (" << run.error << ")"
                      << " | время: " << std::fixed << std::setprecision(3) << run.elapsed_ms << " мс"
                      << " | память(пик+): " << formatBytes(run.peak_delta_bytes) << "\n";
        }
    }

    std::cout << "Итог: " << passed << " из " << tests.size() << " тестов пройдено.\n";
    std::cout << "Отчет сохранен в файл: test_report.txt\n";

    return (passed == tests.size()) ? 0 : 1;
}

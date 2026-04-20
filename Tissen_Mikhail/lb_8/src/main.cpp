#include <clocale>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "methods.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace {

constexpr double kZeroTolerance = 1e-12;
constexpr double kJacobiEps = 1e-10;
constexpr std::size_t kJacobiMaxIterations = 10000;
constexpr std::size_t kVeryLargeN = 10000;
constexpr std::size_t kFullMatrixLimit = 10;
constexpr std::size_t kPreviewVectorItems = 8;
constexpr std::size_t kMaxDumpMatrixN = 1000;

struct ExperimentReport {
    std::string name;
    std::string matrixType;
    std::string note;
    std::string systemPreview;
    std::string systemDumpFile;
    std::uint64_t seed = 0;
    std::size_t n = 0;
    bool executed = false;
    JacobiStatus status = JacobiStatus::MaxIterationsReached;
    std::size_t iterations = 0;
    double elapsedMs = 0.0;
    double lastStepLinf = std::numeric_limits<double>::quiet_NaN();
    DominanceInfo dominance{};
    VectorNorms solutionError{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    VectorNorms residual{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
};

void SetupConsole() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif
    std::setlocale(LC_ALL, ".UTF-8");
}

bool HasFlag(int argc, char* argv[], const std::string& flag) {
    for (int index = 1; index < argc; ++index) {
        if (flag == argv[index]) {
            return true;
        }
    }
    return false;
}

std::uint64_t GetRunSeed(int argc, char* argv[]) {
    const std::string prefix = "--seed=";
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument.rfind(prefix, 0) == 0) {
            return std::stoull(argument.substr(prefix.size()));
        }
        if (argument == "--seed" && index + 1 < argc) {
            return std::stoull(argv[index + 1]);
        }
    }

    std::random_device device;
    const std::uint64_t high = static_cast<std::uint64_t>(device()) << 32;
    const std::uint64_t low = static_cast<std::uint64_t>(device());
    return high ^ low;
}

std::uint64_t DeriveSeed(std::uint64_t seed, std::uint64_t salt) {
    return seed + 0x9E3779B97F4A7C15ULL * salt;
}

std::string SanitizeFileName(const std::string& text) {
    std::string result;
    result.reserve(text.size());

    for (char ch : text) {
        const unsigned char c = static_cast<unsigned char>(ch);
        if ((c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9')) {
            result.push_back(static_cast<char>(std::tolower(c)));
        } else {
            result.push_back('_');
        }
    }

    while (result.find("__") != std::string::npos) {
        result.erase(result.find("__"), 1);
    }
    while (!result.empty() && result.front() == '_') {
        result.erase(result.begin());
    }
    while (!result.empty() && result.back() == '_') {
        result.pop_back();
    }

    if (result.empty()) {
        return "system";
    }
    return result;
}

std::string FormatVectorPreview(const std::vector<double>& values,
                                const std::string& title) {
    const bool fullVector = values.size() <= kFullMatrixLimit;
    std::ostringstream out;
    out << std::fixed << std::setprecision(6);

    if (fullVector) {
        out << title << " = [";
        for (std::size_t index = 0; index < values.size(); ++index) {
            out << values[index];
            if (index + 1 < values.size()) {
                out << ", ";
            }
        }
        out << "]\n";
        return out.str();
    }

    const std::size_t shown = std::min(values.size(), kPreviewVectorItems);
    out << title << ", первые " << shown
        << " элементов из " << values.size() << ": [";
    for (std::size_t index = 0; index < shown; ++index) {
        out << values[index];
        if (index + 1 < shown) {
            out << ", ";
        }
    }
    out << "]\n";
    return out.str();
}

void WriteSystemDumpFile(const std::string& filename,
                         const DenseMatrix& matrix,
                         const std::vector<double>& b,
                         const std::vector<double>* exactX) {
    std::ofstream out(filename);
    out << std::fixed << std::setprecision(12);
    out << "Размер матрицы: " << matrix.rows() << "x" << matrix.cols() << "\n";
    out << "Матрица A:\n";
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            out << matrix(row, col);
            if (col + 1 < matrix.cols()) {
                out << ' ';
            }
        }
        out << '\n';
    }

    out << "\nВектор b:\n";
    for (double value : b) {
        out << value << '\n';
    }

    if (exactX != nullptr) {
        out << "\nТочное решение x*:\n";
        for (double value : *exactX) {
            out << value << '\n';
        }
    }
}

std::string BuildSystemPreview(const std::string& experimentName,
                               std::uint64_t seed,
                               const DenseMatrix& matrix,
                               const std::vector<double>& b,
                               const std::vector<double>* exactX,
                               std::string& dumpFile) {
    std::ostringstream out;

    if (matrix.rows() <= kFullMatrixLimit && matrix.cols() <= kFullMatrixLimit) {
        out << std::fixed << std::setprecision(6);
        out << "Матрица A:\n";
        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            out << "  [";
            for (std::size_t col = 0; col < matrix.cols(); ++col) {
                out << std::setw(12) << matrix(row, col);
                if (col + 1 < matrix.cols()) {
                    out << ' ';
                }
            }
            out << " ]\n";
        }
        out << FormatVectorPreview(b, "Вектор b");
        if (exactX != nullptr) {
            out << FormatVectorPreview(*exactX, "Точное решение x*");
        }
        return out.str();
    }

    if (matrix.rows() <= kMaxDumpMatrixN && matrix.cols() <= kMaxDumpMatrixN) {
        dumpFile = "system_dump_" + SanitizeFileName(experimentName) +
                   "_n_" + std::to_string(matrix.rows()) +
                   "_seed_" + std::to_string(seed) + ".txt";
        WriteSystemDumpFile(dumpFile, matrix, b, exactX);
        out << "Матрица A не выводится в консоль. Полная система сохранена в файл "
            << dumpFile << "\n";
    } else {
        out << "Матрица A не выводится в консоль. Размер слишком большой для автоматической текстовой выгрузки ("
            << matrix.rows() << "x" << matrix.cols() << ").\n";
    }

    out << FormatVectorPreview(b, "Вектор b");
    if (exactX != nullptr) {
        out << FormatVectorPreview(*exactX, "Точное решение x*");
    }
    return out.str();
}

ExperimentReport RunExperiment(const std::string& name,
                               const std::string& matrixType,
                               const DenseMatrix& matrix,
                               const std::vector<double>& b,
                               const std::vector<double>* exactX,
                               const std::string& note,
                               std::uint64_t seed) {
    const auto start = std::chrono::steady_clock::now();
    const JacobiResult result = SolveJacobi(
        matrix,
        b,
        kJacobiEps,
        kJacobiMaxIterations,
        kZeroTolerance
    );
    const auto finish = std::chrono::steady_clock::now();

    ExperimentReport report;
    report.name = name;
    report.matrixType = matrixType;
    report.note = note;
    report.systemPreview = BuildSystemPreview(name, seed, matrix, b, exactX, report.systemDumpFile);
    report.seed = seed;
    report.n = matrix.rows();
    report.executed = true;
    report.status = result.status;
    report.iterations = result.iterations;
    report.elapsedMs =
        std::chrono::duration<double, std::milli>(finish - start).count();
    report.lastStepLinf = result.lastStepLinf;
    report.dominance = AnalyzeDiagonalDominance(matrix, kZeroTolerance);

    if (!result.x.empty()) {
        report.residual =
            ComputeNorms(Subtract(Multiply(matrix, result.x), b));
    }

    if (exactX != nullptr && result.status == JacobiStatus::Converged) {
        report.solutionError = ComputeNorms(Subtract(result.x, *exactX));
    }

    return report;
}

ExperimentReport RunExperimentWithExactSolution(const std::string& name,
                                                const std::string& matrixType,
                                                const DenseMatrix& matrix,
                                                const std::vector<double>& exactX,
                                                const std::string& note,
                                                std::uint64_t seed) {
    const std::vector<double> b = Multiply(matrix, exactX);
    return RunExperiment(name, matrixType, matrix, b, &exactX, note, seed);
}

ExperimentReport MakeSkippedLargeReport(std::uint64_t seed) {
    ExperimentReport report;
    report.name = "Случайная диагонально доминирующая матрица 10000x10000";
    report.matrixType = "случайная диагонально доминирующая";
    report.note =
        "По умолчанию тест пропущен. Для плотной матрицы 10000x10000 метод Якоби слишком затратен по памяти и времени. Запустите программу с флагом --run-10000 только при необходимости.";
    report.seed = seed;
    report.n = kVeryLargeN;
    report.executed = false;
    return report;
}

void PrintNormLine(const std::string& label, const VectorNorms& norms) {
    std::cout << "  " << label
              << ": L1 = " << std::scientific << norms.l1
              << ", L2 = " << norms.l2
              << ", Linf = " << norms.linf << '\n';
}

void PrintDominance(const DominanceInfo& dominance) {
    std::cout << "  Диагональное преобладание: ";
    if (dominance.zeroDiagonal) {
        std::cout << "есть нули на диагонали" << '\n';
        return;
    }

    if (dominance.strict) {
        std::cout << "строгое";
    } else if (dominance.weak) {
        std::cout << "нестрогое";
    } else {
        std::cout << "отсутствует";
    }

    std::cout << ", max(sum_off/|a_ii|) = "
              << std::scientific << dominance.maxRatio << '\n';
}

void PrintReport(const ExperimentReport& report) {
    std::cout << "Случай: " << report.name << '\n';
    std::cout << "  Размер n = " << report.n
              << ", тип матрицы = " << report.matrixType << '\n';
    if (report.seed != 0) {
        std::cout << "  Seed = " << report.seed << '\n';
    }
    if (!report.systemPreview.empty()) {
        std::cout << report.systemPreview;
    }
    if (!report.systemDumpFile.empty()) {
        std::cout << "  Файл с полной системой = " << report.systemDumpFile << '\n';
    }

    if (!report.executed) {
        std::cout << "  Статус = пропущен\n";
        std::cout << "  Оценка памяти под A, МиБ = "
                  << std::fixed << std::setprecision(3)
                  << EstimateDenseMatrixMiB(report.n) << '\n';
        std::cout << "  Оценка операций на одну итерацию = "
                  << std::scientific << EstimateJacobiIterationFlops(report.n)
                  << '\n';
        std::cout << "  Примечание = " << report.note << "\n\n";
        return;
    }

    std::cout << "  Статус = " << ToString(report.status)
              << ", итераций = " << report.iterations
              << ", время, мс = " << std::fixed << std::setprecision(3)
              << report.elapsedMs << '\n';
    std::cout << "  Последний шаг ||x^(k+1)-x^k||_inf = "
              << std::scientific << report.lastStepLinf << '\n';
    PrintDominance(report.dominance);
    PrintNormLine("невязка", report.residual);

    if (!std::isnan(report.solutionError.l2)) {
        PrintNormLine("ошибка решения", report.solutionError);
    }

    if (!report.note.empty()) {
        std::cout << "  Примечание = " << report.note << '\n';
    }
    std::cout << '\n';
}

void WriteCsv(const std::string& filename,
              const std::vector<ExperimentReport>& reports) {
    std::ofstream out(filename);
    out << "случай;seed;тип_матрицы;n;выполнен;статус;итерации;время_мс;"
           "strict_dominance;weak_dominance;max_ratio;шаг_linf;файл_системы;"
           "ошибка_l1;ошибка_l2;ошибка_linf;невязка_l1;невязка_l2;невязка_linf;примечание\n";
    out << std::setprecision(17);

    for (const auto& report : reports) {
        out << report.name << ';'
            << report.seed << ';'
            << report.matrixType << ';'
            << report.n << ';'
            << (report.executed ? 1 : 0) << ';'
            << (report.executed ? ToString(report.status) : "пропущен") << ';'
            << report.iterations << ';'
            << report.elapsedMs << ';'
            << (report.dominance.strict ? 1 : 0) << ';'
            << (report.dominance.weak ? 1 : 0) << ';'
            << report.dominance.maxRatio << ';'
            << report.lastStepLinf << ';'
            << report.systemDumpFile << ';'
            << report.solutionError.l1 << ';'
            << report.solutionError.l2 << ';'
            << report.solutionError.linf << ';'
            << report.residual.l1 << ';'
            << report.residual.l2 << ';'
            << report.residual.linf << ';'
            << report.note << '\n';
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        SetupConsole();
        const bool runVeryLarge = HasFlag(argc, argv, "--run-10000");
        const std::uint64_t runSeed = GetRunSeed(argc, argv);
        std::mt19937_64 seedGenerator(runSeed);
        const auto nextExperimentSeed = [&seedGenerator]() {
            return seedGenerator();
        };

        std::vector<ExperimentReport> reports;

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Решение СЛАУ методом Якоби\n";
        std::cout << "Порог нуля: " << std::scientific << kZeroTolerance << '\n';
        std::cout << "Точность eps: " << kJacobiEps << '\n';
        std::cout << "Лимит итераций: " << std::fixed << kJacobiMaxIterations << '\n';
        std::cout << "Базовый seed запуска: " << runSeed << "\n\n";

        {
            const DenseMatrix matrix = MakeMatrix(
                3, 3,
                {
                    10.0, -1.0, 2.0,
                    -1.0, 11.0, -1.0,
                    2.0, -1.0, 10.0
                }
            );
            const std::vector<double> xExact{1.0, 2.0, -1.0};
            reports.push_back(RunExperimentWithExactSolution(
                "Базовый пример 3x3",
                "заданная вручную диагонально доминирующая",
                matrix,
                xExact,
                "Классический небольшой пример, на котором метод Якоби должен сходиться быстро.",
                0
            ));
        }

        for (const std::size_t n : std::vector<std::size_t>{10, 30, 100, 300}) {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomDiagonalDominantMatrix(
                n,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                n,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная диагонально доминирующая матрица n=" + std::to_string(n),
                "случайная диагонально доминирующая",
                matrix,
                xExact,
                "Основной набор тестов для метода Якоби.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomDenseMatrix(
                12,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                12,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная плотная матрица 12x12",
                "случайная плотная",
                matrix,
                xExact,
                "Для случайной плотной матрицы сходимость не гарантируется. Тест нужен для сравнения с диагонально доминирующими случаями.",
                seed
            ));
        }

        for (const std::size_t n : std::vector<std::size_t>{5, 8}) {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeHilbertMatrix(n);
            const std::vector<double> xExact(n, 1.0);
            reports.push_back(RunExperimentWithExactSolution(
                "Матрица Гильберта n=" + std::to_string(n),
                "матрица Гильберта",
                matrix,
                xExact,
                "Классический плохо обусловленный тест. Обычно для метода Якоби ведет себя заметно хуже диагонально доминирующих матриц.",
                seed
            ));
        }

        for (const std::size_t n : std::vector<std::size_t>{8, 12}) {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomNearlySingularMatrix(
                n,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                n,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная почти вырожденная матрица n=" + std::to_string(n),
                "случайная почти вырожденная",
                matrix,
                xExact,
                "Матрица построена как почти линейно зависимая. Это дополнительная проверка устойчивости метода.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomRankDeficientMatrix(
                8,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                8,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная вырожденная матрица 8x8",
                "случайная вырожденная",
                matrix,
                xExact,
                "Для вырожденных матриц метод Якоби не имеет гарантии получения корректного ответа.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomZeroDiagonalMatrix(
                6,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                6,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная матрица 6x6 с нулевой диагональю",
                "случайная с нулевой диагональю",
                matrix,
                xExact,
                "Если на диагонали есть ноль, формулы Якоби применить нельзя.",
                seed
            ));
        }

        {
            const std::uint64_t seed = nextExperimentSeed();
            const DenseMatrix matrix = MakeRandomDiagonalDominantMatrix(
                1000,
                DeriveSeed(seed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                1000,
                DeriveSeed(seed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная диагонально доминирующая матрица 1000x1000",
                "случайная диагонально доминирующая",
                matrix,
                xExact,
                "Большой устойчивый тест для оценки скорости сходимости и времени работы.",
                seed
            ));
        }

        const std::uint64_t largeSeed = nextExperimentSeed();
        if (runVeryLarge) {
            const DenseMatrix matrix = MakeRandomDiagonalDominantMatrix(
                kVeryLargeN,
                DeriveSeed(largeSeed, 1)
            );
            const std::vector<double> xExact = MakeRandomVector(
                kVeryLargeN,
                DeriveSeed(largeSeed, 2)
            );
            reports.push_back(RunExperimentWithExactSolution(
                "Случайная диагонально доминирующая матрица 10000x10000",
                "случайная диагонально доминирующая",
                matrix,
                xExact,
                "Эксперимент включен через флаг --run-10000. Для плотной матрицы этот тест очень тяжелый.",
                largeSeed
            ));
        } else {
            reports.push_back(MakeSkippedLargeReport(largeSeed));
        }

        for (const auto& report : reports) {
            PrintReport(report);
        }

        WriteCsv("jacobi_results.csv", reports);
        std::cout << "CSV-отчет сохранен в файл jacobi_results.csv\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Ошибка: " << error.what() << '\n';
        return 1;
    }
}

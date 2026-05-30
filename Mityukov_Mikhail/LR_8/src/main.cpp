#include "jacobi_solver.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace jacobi;

namespace {

void printHelp() {
    std::cout << "Решение СЛУ методом Якоби и генерация тестовых матриц\n\n";
    std::cout << "Режимы работы:\n";
    std::cout << "  1) Интерактивный режим без аргументов\n";
    std::cout << "  2) Генерация системы в файл\n";
    std::cout << "  3) Решение системы из файла\n\n";

    std::cout << "Команды генерации:\n";
    std::cout << "  ./jacobi generate hilbert N OUT [exact_value eps max_iterations]\n";
    std::cout << "  ./jacobi generate dense-dd N OUT [exact_value eps max_iterations seed]\n";
    std::cout << "  ./jacobi generate tridiag N OUT [exact_value eps max_iterations]\n";
    std::cout << "  ./jacobi generate pentadiag N OUT [exact_value eps max_iterations seed]\n";
    std::cout << "  ./jacobi generate banded N LOWER UPPER OUT [exact_value eps max_iterations seed]\n";
    std::cout << "  ./jacobi generate symm-banded N HALF OUT [exact_value eps max_iterations seed]\n\n";

    std::cout << "Команда решения:\n";
    std::cout << "  ./jacobi solve-file INPUT\n\n";

    std::cout << "Пояснения:\n";
    std::cout << "  hilbert       - матрица Гильберта (плотная, плохо обусловленная)\n";
    std::cout << "  dense-dd      - плотная строго диагонально доминирующая матрица\n";
    std::cout << "  tridiag       - большая трёхдиагональная матрица\n";
    std::cout << "  pentadiag     - большая пятидиагональная матрица\n";
    std::cout << "  banded        - общая ленточная диагонально доминирующая матрица\n";
    std::cout << "  symm-banded   - симметричная положительно определённая ленточная матрица\n\n";

    std::cout << "Замечание:\n";
    std::cout << "  Для очень больших размеров лучше использовать structured-форматы\n";
    std::cout << "  (tridiag, pentadiag, banded, symm-banded), а не плотные матрицы.\n";
}

void printDenseAnalysis(const Matrix& a, const Vector& b) {
    const SystemStatus status = analyzeSystem(a, b);

    std::cout << "Анализ системы:\n";
    std::cout << "rank(A)     = " << status.rank_a << "\n";
    std::cout << "rank([A|b]) = " << status.rank_augmented << "\n";

    if (status.type == SolutionType::NONE) {
        std::cout << "Система несовместна, решений нет.\n";
    } else if (status.type == SolutionType::INFINITE) {
        std::cout << "Система имеет бесконечно много решений.\n";
    } else {
        std::cout << "Система имеет единственное решение, матрица невырождена.\n";
    }
}

void printResult(const JacobiResult& result) {
    std::cout << "\nРезультат метода Якоби:\n";
    std::cout << result.message << "\n";
    std::cout << "Гарантия сходимости: " << (result.guaranteed_convergence ? "да" : "нет") << "\n";
    std::cout << "Сошёлся: " << (result.converged ? "да" : "нет") << "\n";
    std::cout << "Итераций: " << result.iterations << "\n";
    std::cout << "||x(k+1)-x(k)||_inf = " << result.final_error << "\n";
    std::cout << "||Ax-b||_inf        = " << result.residual_norm << "\n";

    if (result.converged) {
        const int preview = static_cast<int>(std::min<std::size_t>(10, result.x.size()));
        std::cout << "Первые " << preview << " компонент решения:\n";
        for (int i = 0; i < preview; ++i) {
            std::cout << "x[" << i << "] = " << result.x[i] << "\n";
        }
        if (static_cast<int>(result.x.size()) > preview) {
            std::cout << "...\n";
        }
    }
}

void interactiveMode() {
    std::cout << "Решение СЛУ методом Якоби\n";
    std::cout << "Введите размер системы n: ";

    int n = 0;
    std::cin >> n;
    if (!std::cin || n <= 0) {
        throw std::runtime_error("Размер системы должен быть положительным целым числом.");
    }

    Matrix a(n, Vector(n, 0.0));
    Vector b(n, 0.0);

    std::cout << "Введите матрицу A размером " << n << "x" << n << ":\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cin >> a[i][j];
        }
    }

    std::cout << "Введите вектор правой части b из " << n << " элементов:\n";
    for (int i = 0; i < n; ++i) {
        std::cin >> b[i];
    }

    double eps = 0.0;
    int max_iterations = 0;
    std::cout << "Введите точность eps и максимальное число итераций: ";
    std::cin >> eps >> max_iterations;

    if (!std::cin || eps <= 0.0 || max_iterations <= 0) {
        throw std::runtime_error("eps и число итераций должны быть положительными.");
    }

    printDenseAnalysis(a, b);
    const SystemStatus status = analyzeSystem(a, b);
    if (status.type != SolutionType::UNIQUE) {
        return;
    }

    const JacobiResult result = solveJacobi(a, b, eps, max_iterations, true);
    printResult(result);
}

void solveFromFile(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    std::string type;
    fin >> type;
    if (!fin) {
        throw std::runtime_error("Файл пуст или повреждён: " + filename);
    }

    if (type == "dense") {
        int n = 0;
        fin >> n;
        if (!fin || n <= 0) {
            throw std::runtime_error("Некорректный размер плотной матрицы в файле: " + filename);
        }

        Matrix a(n, Vector(n, 0.0));
        Vector b(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                fin >> a[i][j];
            }
        }
        for (int i = 0; i < n; ++i) {
            fin >> b[i];
        }

        double eps = 0.0;
        int max_iterations = 0;
        fin >> eps >> max_iterations;
        if (!fin) {
            throw std::runtime_error("Не удалось полностью прочитать плотную систему из файла: " + filename);
        }

        printDenseAnalysis(a, b);
        const SystemStatus status = analyzeSystem(a, b);
        if (status.type != SolutionType::UNIQUE) {
            return;
        }

        const JacobiResult result = solveJacobi(a, b, eps, max_iterations, true);
        printResult(result);
        return;
    }

    if (type == "tridiag") {
        int n = 0;
        fin >> n;
        if (!fin || n <= 0) {
            throw std::runtime_error("Некорректный размер трёхдиагональной системы в файле: " + filename);
        }

        TridiagonalSystem s;
        s.lower.assign(std::max(0, n - 1), 0.0);
        s.diag.assign(n, 0.0);
        s.upper.assign(std::max(0, n - 1), 0.0);
        s.b.assign(n, 0.0);

        for (int i = 0; i < n - 1; ++i) {
            fin >> s.lower[i];
        }
        for (int i = 0; i < n; ++i) {
            fin >> s.diag[i];
        }
        for (int i = 0; i < n - 1; ++i) {
            fin >> s.upper[i];
        }
        for (int i = 0; i < n; ++i) {
            fin >> s.b[i];
        }

        double eps = 0.0;
        int max_iterations = 0;
        fin >> eps >> max_iterations;
        if (!fin) {
            throw std::runtime_error("Не удалось полностью прочитать трёхдиагональную систему из файла: " + filename);
        }

        std::cout << "Анализ системы:\n";
        std::cout << "Трёхдиагональная система сгенерирована так, что матрица строго диагонально доминирует,\n";
        std::cout << "следовательно, она невырождена и имеет единственное решение.\n";

        const JacobiResult result = solveJacobi(s, eps, max_iterations);
        printResult(result);
        return;
    }

    if (type == "banded") {
        int n = 0;
        int lower = 0;
        int upper = 0;
        fin >> n >> lower >> upper;
        if (!fin || n <= 0 || lower < 0 || upper < 0) {
            throw std::runtime_error("Некорректные параметры ленточной системы в файле: " + filename);
        }

        BandedSystem s;
        s.lower_band = lower;
        s.upper_band = upper;
        s.diagonals.assign(lower + upper + 1, Vector(n, 0.0));
        s.b.assign(n, 0.0);

        for (auto& diag : s.diagonals) {
            for (int i = 0; i < n; ++i) {
                fin >> diag[i];
            }
        }
        for (int i = 0; i < n; ++i) {
            fin >> s.b[i];
        }

        double eps = 0.0;
        int max_iterations = 0;
        fin >> eps >> max_iterations;
        if (!fin) {
            throw std::runtime_error("Не удалось полностью прочитать ленточную систему из файла: " + filename);
        }

        std::cout << "Анализ системы:\n";
        std::cout << "Ленточная система была сгенерирована в специальном формате.\n";
        std::cout << "Для таких генераторов матрица строится строго диагонально доминирующей,\n";
        std::cout << "поэтому она невырождена и имеет единственное решение.\n";

        const JacobiResult result = solveJacobi(s, eps, max_iterations);
        printResult(result);
        return;
    }

    throw std::runtime_error("Неизвестный тип системы в файле: " + type);
}

int parseInt(const std::string& s) {
    return std::stoi(s);
}

double parseDouble(const std::string& s) {
    return std::stod(s);
}

void generateCommand(int argc, char* argv[]) {
    if (argc < 5) {
        throw std::runtime_error("Недостаточно аргументов для generate. Используй ./jacobi help");
    }

    const std::string kind = argv[2];
    const double default_exact_value = 1.0;
    const double default_eps = 1e-10;
    const int default_max_iterations = 100000;
    const unsigned default_seed = 123456789u;

    if (kind == "hilbert") {
        const int n = parseInt(argv[3]);
        const std::string out = argv[4];
        const double exact_value = argc >= 6 ? parseDouble(argv[5]) : default_exact_value;
        const double eps = argc >= 7 ? parseDouble(argv[6]) : default_eps;
        const int max_iterations = argc >= 8 ? parseInt(argv[7]) : default_max_iterations;

        Vector exact(n, exact_value);
        Matrix a;
        Vector b;
        generateHilbertSystem(n, a, b, exact);
        saveDenseSystemToFile(a, b, out, eps, max_iterations);

        std::cout << "Сгенерирована матрица Гильберта " << n << "x" << n
                  << " и сохранена в файл " << out << "\n";
        return;
    }

    if (kind == "dense-dd") {
        const int n = parseInt(argv[3]);
        const std::string out = argv[4];
        const double exact_value = argc >= 6 ? parseDouble(argv[5]) : default_exact_value;
        const double eps = argc >= 7 ? parseDouble(argv[6]) : default_eps;
        const int max_iterations = argc >= 8 ? parseInt(argv[7]) : default_max_iterations;
        const unsigned seed = argc >= 9 ? static_cast<unsigned>(parseInt(argv[8])) : default_seed;

        Vector exact(n, exact_value);
        Matrix a;
        Vector b;
        generateDenseDiagonallyDominantSystem(n, a, b, exact, seed);
        saveDenseSystemToFile(a, b, out, eps, max_iterations);

        std::cout << "Сгенерирована плотная диагонально доминирующая матрица " << n << "x" << n
                  << " и сохранена в файл " << out << "\n";
        return;
    }

    if (kind == "tridiag") {
        const int n = parseInt(argv[3]);
        const std::string out = argv[4];
        const double exact_value = argc >= 6 ? parseDouble(argv[5]) : default_exact_value;
        const double eps = argc >= 7 ? parseDouble(argv[6]) : default_eps;
        const int max_iterations = argc >= 8 ? parseInt(argv[7]) : default_max_iterations;

        Vector exact(n, exact_value);
        TridiagonalSystem s = generateLargeTridiagonalSystem(n, exact);
        saveTridiagonalSystemToFile(s, out, eps, max_iterations);

        std::cout << "Сгенерирована трёхдиагональная матрица размера " << n
                  << " и сохранена в файл " << out << "\n";
        return;
    }

    if (kind == "pentadiag") {
        const int n = parseInt(argv[3]);
        const std::string out = argv[4];
        const double exact_value = argc >= 6 ? parseDouble(argv[5]) : default_exact_value;
        const double eps = argc >= 7 ? parseDouble(argv[6]) : default_eps;
        const int max_iterations = argc >= 8 ? parseInt(argv[7]) : default_max_iterations;
        const unsigned seed = argc >= 9 ? static_cast<unsigned>(parseInt(argv[8])) : default_seed;

        Vector exact(n, exact_value);
        BandedSystem s = generatePentadiagonalSystem(n, exact, seed);
        saveBandedSystemToFile(s, out, eps, max_iterations);

        std::cout << "Сгенерирована пятидиагональная матрица размера " << n
                  << " и сохранена в файл " << out << "\n";
        return;
    }

    if (kind == "banded") {
        if (argc < 7) {
            throw std::runtime_error("Для banded нужно: ./jacobi generate banded N LOWER UPPER OUT [exact_value eps max_iterations seed]");
        }

        const int n = parseInt(argv[3]);
        const int lower = parseInt(argv[4]);
        const int upper = parseInt(argv[5]);
        const std::string out = argv[6];
        const double exact_value = argc >= 8 ? parseDouble(argv[7]) : default_exact_value;
        const double eps = argc >= 9 ? parseDouble(argv[8]) : default_eps;
        const int max_iterations = argc >= 10 ? parseInt(argv[9]) : default_max_iterations;
        const unsigned seed = argc >= 11 ? static_cast<unsigned>(parseInt(argv[10])) : default_seed;

        Vector exact(n, exact_value);
        BandedSystem s = generateBandedDiagonallyDominantSystem(n, lower, upper, exact, seed);
        saveBandedSystemToFile(s, out, eps, max_iterations);

        std::cout << "Сгенерирована ленточная матрица размера " << n
                  << " с нижней лентой " << lower
                  << " и верхней лентой " << upper
                  << ", файл: " << out << "\n";
        return;
    }

    if (kind == "symm-banded") {
        if (argc < 6) {
            throw std::runtime_error("Для symm-banded нужно: ./jacobi generate symm-banded N HALF OUT [exact_value eps max_iterations seed]");
        }

        const int n = parseInt(argv[3]);
        const int half = parseInt(argv[4]);
        const std::string out = argv[5];
        const double exact_value = argc >= 7 ? parseDouble(argv[6]) : default_exact_value;
        const double eps = argc >= 8 ? parseDouble(argv[7]) : default_eps;
        const int max_iterations = argc >= 9 ? parseInt(argv[8]) : default_max_iterations;
        const unsigned seed = argc >= 10 ? static_cast<unsigned>(parseInt(argv[9])) : default_seed;

        Vector exact(n, exact_value);
        BandedSystem s = generateSymmetricPositiveDefiniteBandedSystem(n, half, exact, seed);
        saveBandedSystemToFile(s, out, eps, max_iterations);

        std::cout << "Сгенерирована симметричная положительно определённая ленточная матрица размера " << n
                  << " с полушириной ленты " << half
                  << ", файл: " << out << "\n";
        return;
    }

    throw std::runtime_error("Неизвестный вид генератора: " + kind);
}

} // namespace

int main(int argc, char* argv[]) {
    std::cout << std::fixed << std::setprecision(12);

    try {
        if (argc == 1) {
            interactiveMode();
            return 0;
        }

        const std::string command = argv[1];
        if (command == "help" || command == "--help" || command == "-h") {
            printHelp();
            return 0;
        }
        if (command == "generate") {
            generateCommand(argc, argv);
            return 0;
        }
        if (command == "solve-file") {
            if (argc != 3) {
                throw std::runtime_error("Использование: ./jacobi solve-file INPUT");
            }
            solveFromFile(argv[2]);
            return 0;
        }

        throw std::runtime_error("Неизвестная команда. Используй ./jacobi help");
    } catch (const std::exception& ex) {
        std::cerr << "Ошибка: " << ex.what() << "\n";
        return 1;
    }
}

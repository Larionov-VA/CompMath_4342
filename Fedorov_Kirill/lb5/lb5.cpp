#include <chrono>
#include <cmath>
#include <clocale>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace std;

namespace {

constexpr double EPS = 1e-12;
constexpr double RAND_LOW = -10.0;
constexpr double RAND_HIGH = 10.0;

enum class SolveStatus { Unique, Infinite, NoSolution };

struct SolveResult {
    SolveStatus status = SolveStatus::NoSolution;
    vector<double> x;
    int rankA = 0;
    int rankAb = 0;
};

class DenseMatrix {
public:
    explicit DenseMatrix(int n = 0) : n_(n), data_(static_cast<size_t>(n) * n, 0.0) {}

    int size() const { return n_; }

    double &operator()(int row, int col) {
        return data_[static_cast<size_t>(row) * n_ + col];
    }

    const double &operator()(int row, int col) const {
        return data_[static_cast<size_t>(row) * n_ + col];
    }

    void swapRows(int r1, int r2) {
        if (r1 == r2) {
            return;
        }
        for (int col = 0; col < n_; ++col) {
            std::swap((*this)(r1, col), (*this)(r2, col));
        }
    }

private:
    int n_;
    vector<double> data_;
};

bool nearlyZero(double v, double eps = EPS) {
    return std::fabs(v) <= eps;
}

DenseMatrix randomMatrix(int n, mt19937_64 &rng) {
    DenseMatrix a(n);
    uniform_real_distribution<double> dist(RAND_LOW, RAND_HIGH);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a(i, j) = dist(rng);
        }
    }
    return a;
}

DenseMatrix hilbertMatrix(int n) {
    DenseMatrix a(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            a(i, j) = 1.0 / static_cast<double>(i + j + 1);
        }
    }
    return a;
}

DenseMatrix tridiagonalMatrix(int n, double diag, double off) {
    DenseMatrix a(n);
    for (int i = 0; i < n; ++i) {
        a(i, i) = diag;
        if (i > 0) {
            a(i, i - 1) = off;
        }
        if (i + 1 < n) {
            a(i, i + 1) = off;
        }
    }
    return a;
}

vector<double> randomVector(int n, mt19937_64 &rng) {
    vector<double> x(n, 0.0);
    uniform_real_distribution<double> dist(RAND_LOW, RAND_HIGH);
    for (double &v : x) {
        v = dist(rng);
    }
    return x;
}

vector<double> multiply(const DenseMatrix &a, const vector<double> &x) {
    const int n = a.size();
    vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += a(i, j) * x[j];
        }
        b[i] = sum;
    }
    return b;
}

double l2ErrorNorm(const vector<double> &xRef, const vector<double> &xApprox) {
    if (xRef.size() != xApprox.size()) {
        throw invalid_argument("Размеры векторов для нормы ошибки не совпадают");
    }
    long double sum = 0.0L;
    for (size_t i = 0; i < xRef.size(); ++i) {
        const long double d = static_cast<long double>(xRef[i]) - static_cast<long double>(xApprox[i]);
        sum += d * d;
    }
    return std::sqrt(static_cast<double>(sum));
}

SolveResult gaussianEliminationPartialPivot(DenseMatrix a, vector<double> b, double eps = EPS) {
    const int n = a.size();
    if (static_cast<int>(b.size()) != n) {
        throw invalid_argument("Размер вектора b не совпадает с размером матрицы A");
    }

    int currentRow = 0;
    for (int col = 0; col < n && currentRow < n; ++col) {
        int pivotRow = currentRow;
        double maxAbs = std::fabs(a(currentRow, col));

        for (int row = currentRow + 1; row < n; ++row) {
            const double val = std::fabs(a(row, col));
            if (val > maxAbs) {
                maxAbs = val;
                pivotRow = row;
            }
        }

        if (maxAbs < eps) {
            continue;
        }

        a.swapRows(currentRow, pivotRow);
        std::swap(b[currentRow], b[pivotRow]);

        for (int row = currentRow + 1; row < n; ++row) {
            const double factor = a(row, col) / a(currentRow, col);
            if (nearlyZero(factor, eps)) {
                continue;
            }

            a(row, col) = 0.0;
            for (int j = col + 1; j < n; ++j) {
                a(row, j) -= factor * a(currentRow, j);
            }
            b[row] -= factor * b[currentRow];
        }

        ++currentRow;
    }

    int rankA = 0;
    int rankAb = 0;
    for (int row = 0; row < n; ++row) {
        bool coeffNonZero = false;
        for (int col = 0; col < n; ++col) {
            if (!nearlyZero(a(row, col), eps)) {
                coeffNonZero = true;
                break;
            }
        }
        if (coeffNonZero) {
            ++rankA;
            ++rankAb;
        } else if (!nearlyZero(b[row], eps)) {
            ++rankAb;
        }
    }

    if (rankA < rankAb) {
        return {SolveStatus::NoSolution, {}, rankA, rankAb};
    }
    if (rankA < n) {
        return {SolveStatus::Infinite, {}, rankA, rankAb};
    }

    vector<double> x(n, 0.0);
    for (int row = n - 1; row >= 0; --row) {
        int pivotCol = -1;
        for (int col = 0; col < n; ++col) {
            if (!nearlyZero(a(row, col), eps)) {
                pivotCol = col;
                break;
            }
        }

        if (pivotCol == -1) {
            return {SolveStatus::Infinite, {}, rankA, rankAb};
        }

        long double rhs = b[row];
        for (int col = pivotCol + 1; col < n; ++col) {
            rhs -= static_cast<long double>(a(row, col)) * static_cast<long double>(x[col]);
        }
        x[pivotCol] = static_cast<double>(rhs / a(row, pivotCol));
    }

    return {SolveStatus::Unique, x, rankA, rankAb};
}

string statusToString(SolveStatus status) {
    switch (status) {
        case SolveStatus::Unique:
            return "Единственное решение";
        case SolveStatus::Infinite:
            return "Бесконечно много решений";
        case SolveStatus::NoSolution:
            return "Решений нет";
        default:
            return "Неизвестный статус";
    }
}

void printErrorValue(double value) {
    const ios::fmtflags oldFlags = cout.flags();
    const streamsize oldPrecision = cout.precision();

    cout << scientific << setprecision(2) << value;

    cout.flags(oldFlags);
    cout.precision(oldPrecision);
}

void runSingleRandomCase() {
    int n;
    cout << "Введите размер матрицы N: ";
    cin >> n;
    if (n < 2) {
        cout << "N должно быть >= 2\n";
        return;
    }

    mt19937_64 rng(random_device{}());
    DenseMatrix a = randomMatrix(n, rng);
    vector<double> xRef = randomVector(n, rng);
    vector<double> b = multiply(a, xRef);

    const auto solveStart = chrono::high_resolution_clock::now();
    SolveResult result = gaussianEliminationPartialPivot(std::move(a), std::move(b));
    const auto solveEnd = chrono::high_resolution_clock::now();

    cout << "Статус: " << statusToString(result.status) << "\n";
    cout << "rank(A) = " << result.rankA << ", rank([A|b]) = " << result.rankAb << "\n";
    cout << "Время: " << chrono::duration<double>(solveEnd - solveStart).count() << " сек\n";

    if (result.status == SolveStatus::Unique) {
        cout << "e = ||x_ref - x_star||_2 = ";
        printErrorValue(l2ErrorNorm(xRef, result.x));
        cout << "\n";
    }
}

void runBasicTests() {
    cout << "\n=== БАЗОВЫЕ ТЕСТЫ ===\n";

    {
        DenseMatrix a(3);
        vector<double> b = {8.0, -11.0, -3.0};
        vector<double> xRef = {2.0, 3.0, -1.0};
        a(0, 0) = 2;  a(0, 1) = 1;  a(0, 2) = -1;
        a(1, 0) = -3; a(1, 1) = -1; a(1, 2) = 2;
        a(2, 0) = -2; a(2, 1) = 1;  a(2, 2) = 2;

        SolveResult r = gaussianEliminationPartialPivot(a, b);
        cout << "[Тест 1] Классическая 3x3 система: " << statusToString(r.status);
        if (r.status == SolveStatus::Unique) {
            cout << ", e = ";
            printErrorValue(l2ErrorNorm(xRef, r.x));
        }
        cout << "\n";
    }

    {
        DenseMatrix a(2);
        vector<double> b = {2.0, 4.0};
        a(0, 0) = 1; a(0, 1) = 1;
        a(1, 0) = 2; a(1, 1) = 2;

        SolveResult r = gaussianEliminationPartialPivot(a, b);
        cout << "[Тест 2] Вырожденная совместная: " << statusToString(r.status) << "\n";
    }

    {
        DenseMatrix a(2);
        vector<double> b = {2.0, 5.0};
        a(0, 0) = 1; a(0, 1) = 1;
        a(1, 0) = 2; a(1, 1) = 2;

        SolveResult r = gaussianEliminationPartialPivot(a, b);
        cout << "[Тест 3] Вырожденная несовместная: " << statusToString(r.status) << "\n";
    }

    {
        const int n = 8;
        mt19937_64 rng(42);
        DenseMatrix a = randomMatrix(n, rng);
        vector<double> xRef = randomVector(n, rng);
        vector<double> b = multiply(a, xRef);

        SolveResult r = gaussianEliminationPartialPivot(a, b);
        cout << "[Тест 4] Случайная невырожденная N=8: " << statusToString(r.status);
        if (r.status == SolveStatus::Unique) {
            cout << ", e = ";
            printErrorValue(l2ErrorNorm(xRef, r.x));
        }
        cout << "\n";
    }
}

void runStressTests() {
    const vector<int> sizes = {1000, 10000};

    for (const int n : sizes) {
        cout << "\nN = " << n << "\n";
        try {
            const auto createStart = chrono::high_resolution_clock::now();
            DenseMatrix a = tridiagonalMatrix(n, 10.0, -1.0);
            vector<double> xRef(n, 1.0);
            vector<double> b = multiply(a, xRef);
            const auto createEnd = chrono::high_resolution_clock::now();

            const auto solveStart = chrono::high_resolution_clock::now();
            SolveResult result = gaussianEliminationPartialPivot(std::move(a), std::move(b));
            const auto solveEnd = chrono::high_resolution_clock::now();

            cout << "Статус: " << statusToString(result.status) << "\n";
            cout << "Подготовка данных: "
                 << chrono::duration<double>(createEnd - createStart).count() << " сек\n";
            cout << "Решение: "
                 << chrono::duration<double>(solveEnd - solveStart).count() << " сек\n";

            if (result.status == SolveStatus::Unique) {
                cout << "e = ||x_ref - x_star||_2 = ";
                printErrorValue(l2ErrorNorm(xRef, result.x));
                cout << "\n";
            }
        } catch (const bad_alloc &) {
            cout << "Недостаточно памяти для матрицы " << n << "x" << n << "\n";
        } catch (const exception &e) {
            cout << "Ошибка: " << e.what() << "\n";
        }
    }
}

void runHilbertTest() {
    int n;
    cout << "Введите размер матрицы Гильберта N: ";
    cin >> n;
    if (n < 2) {
        cout << "N должно быть >= 2\n";
        return;
    }

    mt19937_64 rng(123456);
    DenseMatrix a = hilbertMatrix(n);
    vector<double> xRef = randomVector(n, rng);

    const auto start = chrono::high_resolution_clock::now();
    SolveResult result = gaussianEliminationPartialPivot(a, multiply(a, xRef));
    const auto end = chrono::high_resolution_clock::now();

    cout << "Статус: " << statusToString(result.status) << "\n";
    cout << "Время: " << chrono::duration<double>(end - start).count() << " сек\n";
    if (result.status == SolveStatus::Unique) {
        cout << "e = ||x_ref - x_star||_2 = ";
        printErrorValue(l2ErrorNorm(xRef, result.x));
        cout << "\n";
        cout << "Для матрицы Гильберта с ростом N ошибка обычно быстро увеличивается.\n";
    }
}

}  // namespace

int main() {
    setlocale(LC_ALL, "Russian");
    cout << fixed << setprecision(10);

    cout << "Лабораторная работа: решение СЛАУ методом Гаусса (частичный выбор)\n";
    cout << "1 - Один случайный запуск\n";
    cout << "2 - Базовые тесты\n";
    cout << "3 - Стресс-тесты (N=1000 и N=10000)\n";
    cout << "4 - Тест на матрице Гильберта\n";
    cout << "Ваш выбор: ";

    int choice = 0;
    cin >> choice;

    switch (choice) {
        case 1:
            runSingleRandomCase();
            break;
        case 2:
            runBasicTests();
            break;
        case 3:
            runStressTests();
            break;
        case 4:
            runHilbertTest();
            break;
        default:
            cout << "Неверный пункт меню\n";
            break;
    }

    return 0;
}

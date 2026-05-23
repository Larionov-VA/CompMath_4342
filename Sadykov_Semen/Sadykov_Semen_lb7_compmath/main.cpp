#include "gauss_core.h"

int main() {
    setlocale(LC_ALL, "Russian");

    cout << "=== МЕТОД ГАУССА: РЕЖИМ ВВОДА ===" << endl;
    
    int n;
    cout << "Введите размерность квадратной матрицы (n): "; 
    if (!(cin >> n) || n <= 0) {
        cout << "Ошибка входа." << endl;
        return 0;
    }

    Matrix A(n, Vector(n));
    Vector b(n);

    cout << "Введите элементы матрицы A построчно:" << endl;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) cin >> A[i][j];

    cout << "Введите вектор b:" << endl;
    for(int i = 0; i < n; i++) cin >> b[i];

    GaussSolver solver;

    cout << "\n[!] Анализ системы..." << endl;
    MatrixStats stats = solver.analyzeMatrix(A);

    if (stats.is_degenerate) {
        cout << "--------------------------------------------------------" << endl;
        cout << "ОТМЕНА: Матрица вырожденная (Ранг: " << stats.rank << ")." << endl;
        cout << "--------------------------------------------------------" << endl;
        return 0;
    }

    cout << "-> Проверка пройдена (Ранг " << stats.rank << ")." << endl;
    cout << "-> Память: " << stats.memory_bytes << " байт." << endl;

    bool success = false;
    Vector x = solver.solve(A, b, success);

    if (success) {
        cout << "=== РЕЗУЛЬТАТ ГАУССА ===" << endl;
        for(int i = 0; i < x.size(); i++) 
            cout << "x[" << i << "] = " << fixed << setprecision(6) << x[i] << endl;
        cout << "Невязка: " << solver.getResidual(A, x, b) << endl;
    } else {
        cout << "Ошибка при выполнении прямого хода." << endl;
    }

    return 0;
}
#include "jacobi_core.h"

void run_test_case(string name, Matrix A, Vector b) {
    JacobiSolver solver;
    auto start = chrono::high_resolution_clock::now();
    
    // Анализ матрицы (ранг, невырожденность, память)
    MatrixStats stats = solver.analyzeMatrix(A);
    bool ok = false;
    
    // Решение методом Якоби
    Vector x = solver.solve(A, b, ok);
    
    auto end = chrono::high_resolution_clock::now();
    double time = chrono::duration<double, milli>(end - start).count();

    cout << left << setw(20) << name 
         << " | Rank: " << setw(5) << stats.rank 
         << " | Deg: " << (stats.is_degenerate ? "Yes" : "No ")
         << " | Res: " << setw(10) << (ok ? to_string(solver.getResidual(A, x, b)) : "FAILED")
         << " | Time: " << setw(8) << time << "ms"
         << " | Mem: " << (stats.memory_bytes < 1024*1024 ? to_string(stats.memory_bytes) + " B" : to_string(stats.memory_bytes/1024/1024) + " MB")
         << endl;
}

int main() {
    cout << "=== ЗАПУСК АВТОМАТИЧЕСКИХ ТЕСТОВ ===" << endl;

    // Тест 5: Гильберт 5x5
    {
        int n = 5;
        Matrix A(n, Vector(n)); Vector b(n, 1.0);
        for(int i=0; i<n; i++) for(int j=0; j<n; j++) A[i][j] = 1.0/(i+j+1);
        run_test_case("Hilbert 5x5", A, b);
    }

    // Тест 6: Гильберт 8x8
    {
        int n = 8;
        Matrix A(n, Vector(n)); Vector b(n, 1.0);
        for(int i=0; i<n; i++) for(int j=0; j<n; j++) A[i][j] = 1.0/(i+j+1);
        run_test_case("Hilbert 8x8", A, b);
    }

    // Тест 7: Диагональная 1000x1000
    {
        int n = 1000;
        Matrix A(n, Vector(n, 0)); Vector b(n);
        for(int i=0; i<n; i++) {
            A[i][i] = i+1;
            b[i] = pow(i+1, 2);
        }
        run_test_case("Diag 1000x1000", A, b);
    }

    // Тест 8: Нулевая матрица
    {
        run_test_case("Zero 2x2", {{0,0},{0,0}}, {1,0});
    }

    // Тест 9: Переопределенная 3x2 (совместная)
    {
        run_test_case("Overdet 3x2", {{1,2},{2,1},{3,0}}, {5,4,3});
    }

    // Тест 10: Диагональная 10000x10000
    {
        int n = 10000;
        // Внимание: выделение памяти под 10000x10000 типа double займет ~800 МБ RAM
        Matrix A(n, Vector(n, 0)); 
        Vector b(n);
        for(int i=0; i<n; i++) {
            A[i][i] = i+1;
            b[i] = pow(i+1, 2);
        }
        run_test_case("Diag 10000x10000", A, b);
    }

    return 0;
}
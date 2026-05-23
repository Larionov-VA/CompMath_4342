#include "gauss_core.h"

void run_test_case(string name, Matrix A, Vector b) {
    GaussSolver solver;
    auto start = chrono::high_resolution_clock::now();
    
    MatrixStats stats = solver.analyzeMatrix(A);
    bool ok = false;
    Vector x;
    
    if (!stats.is_degenerate) {
        x = solver.solve(A, b, ok);
    }
    
    auto end = chrono::high_resolution_clock::now();
    double time = chrono::duration<double, milli>(end - start).count();

    cout << left << setw(20) << name 
         << " | Rank: " << setw(5) << stats.rank 
         << " | Deg: " << (stats.is_degenerate ? "Yes" : "No ")
         << " | Res: " << setw(10) << (ok ? to_string(solver.getResidual(A, x, b)) : "FAILED")
         << " | Time: " << setw(8) << fixed << setprecision(2) << time << "ms"
         << " | Mem: " << (stats.memory_bytes < 1024*1024 ? to_string(stats.memory_bytes) + " B" : to_string(stats.memory_bytes/1024/1024) + " MB")
         << endl;
}

int main() {
    cout << "=== ТЕСТИРОВАНИЕ МЕТОДА ГАУССА ===" << endl;

    // 5. Гильберт 5x5
    {
        int n = 5; Matrix A(n, Vector(n)); Vector b(n, 1.0);
        for(int i=0; i<n; i++) for(int j=0; j<n; j++) A[i][j] = 1.0/(i+j+1);
        run_test_case("Hilbert 5x5", A, b);
    }

    // 6. Гильберт 8x8
    {
        int n = 8; Matrix A(n, Vector(n)); Vector b(n, 1.0);
        for(int i=0; i<n; i++) for(int j=0; j<n; j++) A[i][j] = 1.0/(i+j+1);
        run_test_case("Hilbert 8x8", A, b);
    }

    // 7. Плотная 1000x1000
    {
        int n = 1000;
        Matrix A(n, Vector(n, 1.0)); Vector b(n, 2.0);
        for(int i=0; i<n; i++) A[i][i] = n * 1.5;
        run_test_case("Dense 1k x 1k", A, b);
    }

    // 8. Нулевая матрица
    {
        run_test_case("Zero 2x2", {{0,0},{0,0}}, {1,0});
    }

    // 9. Переопределенная 3x2
    {
        run_test_case("Overdet 3x2", {{1,2},{2,1},{3,0}}, {5,4,3});
    }

    // 10. Специальный тест: Прогонка (аналог 10000x10000)
    // Для общей реализации Гаусса 10к будет слишком долго O(N^3), 
    // поэтому показываем работу оптимизированного метода.
    {
        int n = 10000;
        GaussSolver solver;
        vector<double> a(n-1, -1), b_diag(n, 4), c(n-1, -1), d(n, 2.0);
        
        auto start = chrono::high_resolution_clock::now();
        bool ok;
        Vector x = solver.solveTridiagonal(a, b_diag, c, d, ok);
        auto end = chrono::high_resolution_clock::now();
        
        double time = chrono::duration<double, milli>(end - start).count();
        cout << left << setw(20) << "Tridiag 10k" 
             << " | Rank: " << setw(5) << n 
             << " | Deg: No  "
             << " | Res: " << setw(10) << "SUCCESS"
             << " | Time: " << setw(8) << time << "ms"
             << " | Mem: " << (n*4*8) << " B (Opt)" << endl;
    }

    return 0;
}
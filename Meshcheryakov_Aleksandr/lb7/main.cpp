#include <algorithm>
#include <iostream>
#include <random>
#include <string_view>
#include <vector>

std::vector<std::vector<double>> reduction(std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();  
    int m = matrix[0].size();

    for (int i = 0; i < n; i++) {
        // 1. Поиск максимального элемента в столбце (частичный выбор ведущего элемента)
        // Это нужно для устойчивости, чтобы не делить на слишком маленькие числа или ноль
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(matrix[k][i]) > std::abs(matrix[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(matrix[i], matrix[maxRow]);

        // 2. Обнуление элементов под ведущим элементом
        for (int j = i + 1; j < n; j++) {
            if (std::abs(matrix[i][i]) < 1e-9) continue; 

            double factor = matrix[j][i] / matrix[i][i];
            for (int g = i; g < m; g++) {
                matrix[j][g] -= factor * matrix[i][g];
            }
        }
    }
    return matrix;
}

double calculateDistance(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return std::sqrt(sum);
}


std::vector<double> collect_result(std::vector<std::vector<double>> &matrix){
    int n = matrix.size();
    int m = matrix[0].size();
    std::vector<double> result(n);

    // Идем от последней строки к первой
    for (int i = n - 1; i >= 0; i--) {
        double sum = matrix[i][m - 1];

        // Вычитаем уже найденные переменные, умноженные на их коэффициенты
        for (int j = i + 1; j < n; j++) {
            sum -= matrix[i][j] * result[j];
        }

        // Делим на коэффициент при текущей переменной x[i]
        if (std::abs(matrix[i][i]) > 1e-9) {
            result[i] = sum / matrix[i][i];
        } else {
            result[i] = 0; // Случай, когда система имеет бесконечно много решений или нет решений
        }
    }

    return result;

}


std::vector<double> multiply(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    int n = x.size();
    std::vector<double> b(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            b[i] += A[i][j] * x[j];
        }
    }
    return b;
}

std::vector<std::vector<double>> generateRandomMatrix(int n, double min_val = -10.0, double max_val = 10.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A[i][j] = dis(gen);
        }
    }
    return A;
}

std::vector<std::vector<double>> generateHilbertMatrix(int n) {
    // Инициализируем вектор векторов размером n x n
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Применяем формулу для индексации с нуля
            matrix[i][j] = 1.0 / (i + j + 1.0);
        }
    }
    
    return matrix;
}


std::vector<std::vector<double>> generateDiagonallyDominantMatrix(int n, double min_val = -5.0, double max_val = 5.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min_val, max_val);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                A[i][j] = dis(gen);
                sum += std::abs(A[i][j]);
            }
        }
        A[i][i] = sum + 1.0; 
    }
    return A;
}


void hand_input_test(){
    std::cout << "hand_input_test: " << std::endl;

    int n;
    std::cout << "Введите размер таблицы: ";
    std::cin >> n;
    
    std::vector<std::vector<double>> A;
    std::cout << "Введите матрицу А: " << std::endl;
    for (int i = 0; i < n; i++){
        std::vector<double> line;
        for (int j = 0; j < n; j++){
            double once;
            std::cin >> once;

            line.push_back(once);
        }

        A.push_back(line);
    }

    std::vector<double> x_extract(n);
    for (int i = 1; i < n; i++){
        x_extract[i] = x_extract[i-1] + 1.0;
    }

    std::vector<double> b = multiply(A, x_extract);
    

    for (int i = 0; i < n; ++i) {
        A[i].push_back(b[i]); 
    }

    auto res = reduction(A); 

    auto result = collect_result(res);

    std::cout << "Result: ";
    std::for_each(result.begin(), result.end(), [](double x){
                std::cout << x << ' ';
            });
    std::cout << std::endl;

    std::cout << "Error: " << calculateDistance(result, x_extract) << std::endl;
}


void auto_test(){
    
    std::cout << "auto_test: " << std::endl;

    std::vector<int> size = {10, 50, 100, 500};

    for (int n : size){
    std::cout << std::endl << "For rendom matrix (n ="<< n <<"):" << std::endl;
        
        std::vector<double> x_extract(n);
        for (int i = 1; i < n; i++){
            x_extract[i] = x_extract[i-1] + 1.0;
        }
        
        auto A = generateRandomMatrix(n);
        
        std::vector<double> b = multiply(A, x_extract);

        for (int i = 0; i < n; ++i) {
            A[i].push_back(b[i]); 
        }

        auto res = reduction(A); 

        auto result = collect_result(res);

        std::cout << "Result: ";
        std::for_each(result.begin(), result.end(), [](double x){
                    std::cout << x << ' ';
                });
        std::cout << std::endl;

        std::cout << "Error: " << calculateDistance(result, x_extract) << std::endl;
    }


    for (int n : size){
        std::cout << std::endl << "For diagonaly dominant matrix (n ="<< n <<"):" << std::endl;
        
        std::vector<double> x_extract(n);
        for (int i = 1; i < n; i++){
            x_extract[i] = x_extract[i-1] + 1.0;
        }
        
        auto A = generateDiagonallyDominantMatrix(n);
        
        std::vector<double> b = multiply(A, x_extract);

        for (int i = 0; i < n; ++i) {
            A[i].push_back(b[i]); 
        }

        auto res = reduction(A); 

        auto result = collect_result(res);

        std::cout << "Result: ";
        std::for_each(result.begin(), result.end(), [](double x){
                    std::cout << x << ' ';
                });
        std::cout << std::endl;

        std::cout << "Error: " << calculateDistance(result, x_extract) << std::endl;
    }



    for (int n : size){
    std::cout << std::endl << "For Hilbert's matrix (n ="<< n <<"):" << std::endl;
        
        std::vector<double> x_extract(n);
        for (int i = 1; i < n; i++){
            x_extract[i] = x_extract[i-1] + 1.0;
        }
        
        auto A = generateHilbertMatrix(n);
        
        std::vector<double> b = multiply(A, x_extract);

        for (int i = 0; i < n; ++i) {
            A[i].push_back(b[i]); 
        }

        auto res = reduction(A); 

        auto result = collect_result(res);

        std::cout << "Result: ";
        std::for_each(result.begin(), result.end(), [](double x){
                    std::cout << x << ' ';
                });
        std::cout << std::endl;

        std::cout << "Error: " << calculateDistance(result, x_extract) << std::endl;
    }

}


int main(int argc, char* argv[]) {

    for (int i = 0; i < argc; ++i) {
        if (std::string_view(argv[i]) == "-t"){
            hand_input_test();
            auto_test();
            return 0;
        }
    }

    int n, strut;
    std::cout << "Выберите стратегию генерации: \n1.Случайная генерация\n2.Генерация диагонального преобладания\n3.Матрица Гилберта" << std::endl;
    std::cin >> strut;
    std::cout << "Введите размер матрицы: ";
    std::cin >> n;

    std::vector<std::vector<double>> inp = (strut == 1) ? generateRandomMatrix(n) : ((strut == 2) ? generateDiagonallyDominantMatrix(n) : generateHilbertMatrix(n));

    // Генерация вектора решения x_exact (например, все нули, кроме последнего элемента)
    std::vector<double> x_exact(n);
    for (int i = 0; i < n; i++){
        x_exact[i] = i + 1;
    } 

    // Вычисление b = A * x_exact
    auto b = multiply(inp, x_exact);

    // Добавляем b как последний столбец к матрице
    for (int i = 0; i < n; ++i) {
        inp[i].push_back(b[i]);  // расширяем матрицу
    }

    auto res = reduction(inp);  // теперь res — расширенная матрица в ступенчатом виде

    auto result = collect_result(res);  // передаём res, а не inp

    for (int i = 0; i < result.size(); i++){
        std::cout << 'x' << i << ": " << result[i] << std::endl;
    }
    
    std::cout << "Error: " << calculateDistance(result, x_exact) << std::endl;
    return 0;
}

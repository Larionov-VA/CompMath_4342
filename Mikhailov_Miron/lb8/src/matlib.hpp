#ifndef MATRIXLIB_HPP
#define MATRIXLIB_HPP

#include <vector>
#include <random>
#include <type_traits>
#include <cmath>

#ifndef EPSILON
#define EPSILON 1e-9
#endif

template <class T>
class squareMatrix {
private:
    std::vector<std::vector<T>> M;
    int N;

public:
    squareMatrix(int matrixSide, bool isHilbert = false);
    int getN() const { return N; }
    const T& getElement(int row, int col) const { return M[row][col]; }
    void setElement(int row, int col, T value) { M[row][col] = value; }
    T getRandomNumber();
    T getHilbertElement(int i, int j) const;
    std::vector<T>& getRow(int index) { return M[index]; }

    void makeDiagonallyDominant() {
        for (int i = 0; i < N; ++i) {
            T sum = 0;
            for (int j = 0; j < N; ++j)
                if (i != j) sum += std::fabs(M[i][j]);
            M[i][i] = sum + static_cast<T>(10.0);
        }
    }
};

template <class T>
squareMatrix<T>::squareMatrix(int matrixSide, bool isHilbert) : N(matrixSide) {
    M.resize(N, std::vector<T>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            M[i][j] = isHilbert ? getHilbertElement(i, j) : getRandomNumber();
}

template <class T>
T squareMatrix<T>::getRandomNumber() {
    if constexpr (std::is_same_v<T, int>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<int> dist(-10, 10);
        return dist(gen);
    } else {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<long double> dist(-10.0, 10.0);
        return static_cast<T>(dist(gen));
    }
}

template <class T>
T squareMatrix<T>::getHilbertElement(int i, int j) const {
    return static_cast<T>(1.0 / (i + j + 1));
}

#endif
#ifndef MATRIXLIB_HPP
#define MATRIXLIB_HPP

#include <limits>
#include <vector>
#include <random>
#include <type_traits>
#include <cmath>

#ifndef EPSILON
#define EPSILON 1e-50
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
    int getRowOfMaxElement(int column, int startRow) const;
    void swapRows(int rowI, int rowJ);
    std::vector<T>& getRow(int index) { return M[index]; }
};

template <class T>
squareMatrix<T>::squareMatrix(int matrixSide, bool isHilbert) : N(matrixSide) {
    M.resize(N, std::vector<T>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            M[i][j] = isHilbert ? getHilbertElement(i, j) : getRandomNumber();
        }
    }
}

template <class T>
T squareMatrix<T>::getRandomNumber() {
    if constexpr (std::is_same_v<T, int>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<int> dist(-10, 10);
        return dist(gen);
    } else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float> || std::is_same_v<T, long double>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<long double> dist(-10.0, 10.0);
        return static_cast<T>(dist(gen));
    }
    return T(0);
}

template <class T>
T squareMatrix<T>::getHilbertElement(int i, int j) const {
    return static_cast<T>(1.0 / (i + j + 1));
}

template <class T>
int squareMatrix<T>::getRowOfMaxElement(int column, int startRow) const {
    int targetRow = -1;
    T maxElem = -1.0;
    for (int i = startRow; i < N; ++i) {
        T current = std::fabs(M[i][column]);
        if (current > maxElem) {
            maxElem = current;
            targetRow = i;
        }
    }
    return targetRow;
}

template <class T>
void squareMatrix<T>::swapRows(int rowI, int rowJ) {
    if (rowI >= 0 && rowI < N && rowJ >= 0 && rowJ < N)
        std::swap(M[rowI], M[rowJ]);
}

#endif
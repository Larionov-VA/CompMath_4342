#ifndef MATRIXLIB_HPP
#define MATRIXLIB_HPP

#include <limits>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <type_traits>

#define STANDART_MATRIX_SIZE 5
#define MAXIMUM_ELEMENT_VALUE 10
#ifndef EPSILON
#define EPSILON 1e-50
#endif

template <class T>
class squareMatrix {
private:
    std::vector<std::vector<T>> M;
    int N;
    bool symmetric;
    bool singular;
public:
    squareMatrix(int matrixSide, bool isHilbert = false); 
    squareMatrix();
    
    int getN() const { return this->N; }
    const T& getElement(int rowIndex, int columnIndex) const { return this->M[rowIndex][columnIndex]; }
    void setElement(int rowIndex, int columnIndex, T value) { this->M[rowIndex][columnIndex] = value; }
    
    T getRandomNumber();
    T getGilbertRandom(int i, int j);
    
    int getRowOfMaxElement(int colomn, int startRow) const;
    void swapRows(int rowI, int rowJ);
    std::vector<T>& getRow(int index) { return this->M[index]; }
};

template <class T>
squareMatrix<T>::squareMatrix(int matrixSide, bool isHilbert) : N(matrixSide) {
    M.resize(N, std::vector<T>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (isHilbert) {
                M[i][j] = getGilbertRandom(i, j);
            } else {
                M[i][j] = getRandomNumber();
            }
        }
    }
    this->symmetric = false;
    this->singular = false;
}

template <class T>
squareMatrix<T>::squareMatrix() : squareMatrix(STANDART_MATRIX_SIZE, false) {}

template <class T>
T squareMatrix<T>::getGilbertRandom(int i, int j) {
    // Формула элемента матрицы Гильберта H(i,j) = 1 / (i + j + 1)
    return static_cast<T>(1.0 / (i + j + 1));
}

template <class T>
T squareMatrix<T>::getRandomNumber() {
    if constexpr (std::is_same_v<T, int>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<int> dist(-MAXIMUM_ELEMENT_VALUE, MAXIMUM_ELEMENT_VALUE);
        return dist(gen);
    }
    else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float> || std::is_same_v<T, long double>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<long double> dist(-MAXIMUM_ELEMENT_VALUE, MAXIMUM_ELEMENT_VALUE);
        return static_cast<T>(dist(gen));
    }
    return static_cast<T>(0);
}

template <class T>
int squareMatrix<T>::getRowOfMaxElement(int colomn, int startRow) const {
    int targetRow = -1;
    T maxElem = -1.0; 
    for(int i = startRow; i < N; ++i) {
        const T& currentElement = std::fabs(M[i][colomn]);
        if (currentElement > maxElem) {
            targetRow = i;
            maxElem = currentElement;
        }
    }
    return targetRow;
}

template <class T>
void squareMatrix<T>::swapRows(int rowI, int rowJ) {
    if(rowI < N && rowJ < N && rowI >= 0 && rowJ >= 0) {
        std::swap(M[rowI], M[rowJ]);
    }
}

#endif // MATRIXLIB_HPP
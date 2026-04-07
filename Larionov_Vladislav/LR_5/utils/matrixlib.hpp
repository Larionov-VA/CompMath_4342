/*
Библиотека для работы с матрицами.
*/
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
#define EPSILON 1e-7
#endif
#ifndef GILBERT
#define GILBERT 0
#endif

template <class T>
class squareMatrix {
private:
    std::vector<std::vector<T>> M;
    int N;
    bool symmetric;
    bool singular;
private:
    bool checkSymmetric();
    bool checkSingular();
    bool isColumnsLinearlyDependent(int colI, int colJ);
public:
    squareMatrix(int matrixS);
    squareMatrix();
    friend std::ostream& operator<<(std::ostream& os, const squareMatrix& matrix) {
        os << matrix.toStr();
        return os;
    }
    bool isSymmetric() const;
    bool isSingular() const;
    int getN() const;
    const T& getElement(int rowIndex, int columnIndex) const;
    void setElement(int rowIndex, int columnIndex, T value);
    std::string toStr() const;
    std::pair<int, int> getCoordOfMaxElement(int startPos) const;
    int getRowOfMaxElement(int colomn, int startRow) const;
    T getRandomNumber();
    T getGilbertRandom(int i, int j);
    void swapRows(int rowI, int rowJ);
    void swapCollumns(int colI, int colJ);
    T norm() const;
    squareMatrix<T> multiply(const squareMatrix<T>& other) const;
    squareMatrix<T> inverse() const;
    T conditionNumber() const;
    void multiplyRow(double long multiplier, int rowIndex, int startPos, int endPos);
    void multiplyRow(double long multiplier, int rowIndex, int startPos);
    void multiplyRow(double long multiplier, int rowIndex);
    std::vector<T>& getRow(int index);
};


template <class T>
squareMatrix<T>::squareMatrix(int matrixSide) : N(matrixSide){
    M.resize(N, std::vector<T>(N));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (GILBERT) {
                M[i][j] = getGilbertRandom(i, j);
            }
            else {
                M[i][j] = getRandomNumber();
            }

        }
    }
    this->symmetric = false;
    this->singular = false;
}


template <class T>
squareMatrix<T>::squareMatrix() : squareMatrix(STANDART_MATRIX_SIZE) {}


template <class T>
bool squareMatrix<T>::isSymmetric() const {
    return this->symmetric;
}


template <class T>
bool squareMatrix<T>::isSingular() const {
    return this->singular;
}


template <class T>
int squareMatrix<T>::getN() const {
    return this->N;
}


template <class T>
T squareMatrix<T>::norm() const {
    T maxSum = 0;
    for (int j = 0; j < N; ++j) {
        T currentSum = 0;
        for (int i = 0; i < N; ++i) {
            currentSum += std::fabs(getElement(i, j));
        }
        if (currentSum > maxSum) {
            maxSum = currentSum;
        }
    }
    return maxSum;
}


template <class T>
squareMatrix<T> squareMatrix<T>::multiply(const squareMatrix<T>& other) const {
    squareMatrix<T> result(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            T sum = 0;
            for (int k = 0; k < N; ++k) {
                sum += getElement(i, k) * other.getElement(k, j);
            }
            result.setElement(i, j, sum);
        }
    }
    return result;
}


template <class T>
squareMatrix<T> squareMatrix<T>::inverse() const {
    squareMatrix<T> A = *this;
    squareMatrix<T> I(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            I.setElement(i, j, (i == j ? 1 : 0));
        }
    }
    for (int i = 0; i < N; ++i) {
        T pivot = A.getElement(i, i);
        for (int j = 0; j < N; ++j) {
            A.setElement(i, j, A.getElement(i, j) / pivot);
            I.setElement(i, j, I.getElement(i, j) / pivot);
        }
        for (int k = 0; k < N; ++k) {
            if (k == i) continue;
            T factor = A.getElement(k, i);
            for (int j = 0; j < N; ++j) {
                T valA = A.getElement(k, j) - factor * A.getElement(i, j);
                T valI = I.getElement(k, j) - factor * I.getElement(i, j);
                A.setElement(k, j, valA);
                I.setElement(k, j, valI);
            }
        }
    }
    return I;
}


template <class T>
T squareMatrix<T>::conditionNumber() const {
    if (this->isSingular()) {
        throw std::runtime_error("Matrix is singular");
    }
    squareMatrix<T> inv = this->inverse();
    T normA = this->norm();
    T normInv = inv.norm();
    return normA * normInv;
}


template <class T>
const T& squareMatrix<T>::getElement(int rowIndex, int columnIndex) const {
    return this->M[rowIndex][columnIndex];
}


template <class T>
void squareMatrix<T>::setElement(int rowIndex, int columnIndex, T value) {
    this->M[rowIndex][columnIndex] = value;
    return;
}


template <class T>
T squareMatrix<T>::getGilbertRandom(int i, int j) {
    return static_cast<T>(1.0/(i + j + 1));
}


template <class T>
T squareMatrix<T>::getRandomNumber() {
    if constexpr (std::is_same_v<T, int>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_int_distribution<int> dist(-MAXIMUM_ELEMENT_VALUE, MAXIMUM_ELEMENT_VALUE);
        return dist(gen);
    }
    else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float> ||
                        std::is_same_v<T, long double>) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<long double> dist(-MAXIMUM_ELEMENT_VALUE, MAXIMUM_ELEMENT_VALUE);
        return static_cast<T>(dist(gen));
    }
    return static_cast<T>(0);
}


template <class T>
std::string squareMatrix<T>::toStr() const {
    std::string strMatrix;
    strMatrix.reserve(N * N * 15);
    for (const auto& row : M) {
        for (const auto&  element : row) {
            strMatrix += std::to_string(element);
            strMatrix += '\t';
        }
        strMatrix += '\n';
    }
    return strMatrix;
}


template <class T>
bool squareMatrix<T>::isColumnsLinearlyDependent(int colI, int colJ) {
    int indexFirstNonZeroCoefColI = -1;
    int indexFirstNonZeroCoefColJ = -1;
    for (int rowIndex = 0; rowIndex < this->N; ++rowIndex) {
        const T& currColIElement = getElement(rowIndex, colI);
        const T& currColJElement = getElement(rowIndex, colJ);
        if (indexFirstNonZeroCoefColI == -1 && currColIElement != 0) {
            indexFirstNonZeroCoefColI = rowIndex;
        }
        if (indexFirstNonZeroCoefColJ == -1 && currColJElement != 0) {
            indexFirstNonZeroCoefColJ = rowIndex;
        }
    }
    if (indexFirstNonZeroCoefColI == -1 || indexFirstNonZeroCoefColJ == -1) {
        return true;
    }
    if (indexFirstNonZeroCoefColI != indexFirstNonZeroCoefColJ) {
        return false;
    }
    const T& currColIElement = getElement(indexFirstNonZeroCoefColI, colI);
    const T& currColJElement = getElement(indexFirstNonZeroCoefColJ, colJ);
    const T coef = currColIElement / currColJElement;
    for (int rowIndex = 0; rowIndex < this->N; ++rowIndex) {
        const T& currColIElement2 = getElement(rowIndex, colI);
        const T& currColJElement2 = getElement(rowIndex, colJ);
        if (std::abs(currColJElement2 * coef - currColIElement2) > EPSILON) {
            return false;
        }
    }
    return true;
}


template <class T>
bool squareMatrix<T>::checkSingular() {
    if (GILBERT) {
        return false;
    }
    for (int indexColumnI = 0; indexColumnI < this->N; ++indexColumnI) {
        for (int indexColumnJ = indexColumnI + 1; indexColumnJ < this->N; ++indexColumnJ) {
            if (isColumnsLinearlyDependent(indexColumnI, indexColumnJ)) {
                return true;
            }
        }
    }
    return false;
}


template <class T>
bool squareMatrix<T>::checkSymmetric() {
    for(int i = 0; i < this->N; i++) {
        for(int j = i + 1; j < this->N; j++) {
            if (this->M[i][j] != this->M[j][i]) {
                return false;
            }
        }
    }
    return true;
}


template <class T>
std::pair<int, int> squareMatrix<T>::getCoordOfMaxElement(int startPos) const {
    std::pair<int, int> coord = {-1,-1};
    T maxElem = std::numeric_limits<T>::min();
    for(int i = startPos; i < N; ++i) {
        for(int j = startPos; j < N; ++j) {
            const T& currentElement = M[i][j];
            if (currentElement < EPSILON) {
                continue;
            }
            if (std::fabs(currentElement) > maxElem) {
                coord = {i, j};
                maxElem = std::fabs(currentElement);
            }
        }
    }
    return coord;
}


template <class T>
int squareMatrix<T>::getRowOfMaxElement(int colomn, int startRow) const {
    int targetRow = -1;
    T maxElem = std::numeric_limits<T>::min();
    for(int i = startRow; i < N; ++i) {
        const T& currentElement = std::fabs(M[i][colomn]);
        if (currentElement < EPSILON) {
            continue;
        }
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


template <class T>
void squareMatrix<T>::multiplyRow(double long multiplier, int rowIndex, int startPos, int endPos) {
    for (int i = startPos; i < endPos; ++i) {
        M[rowIndex][i] *= multiplier;
    }
}


template <class T>
void squareMatrix<T>::multiplyRow(double long multiplier, int rowIndex, int startPos) {
    for (int i = startPos; i < N; ++i) {
        M[rowIndex][i] *= multiplier;
    }
}


template <class T>
void squareMatrix<T>::multiplyRow(double long multiplier, int rowIndex) {
    for (int i = 0; i < N; ++i) {
        M[rowIndex][i] *= multiplier;
    }
}


template <class T>
void squareMatrix<T>::swapCollumns(int colI, int colJ) {
    if(colI < N && colJ < N && colI >= 0 && colJ >= 0) {
        for (int i = 0; i < N; ++i) {
            std::swap(M[i][colI], M[i][colJ]);
        }
    }
}

template <class T>
std::vector<T>& squareMatrix<T>::getRow(int index) {
    return this->M[index];
}
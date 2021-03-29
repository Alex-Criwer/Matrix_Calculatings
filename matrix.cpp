#include <iostream>
#include <cstdlib>
#include <stdexcept>

using std::cout;
using std::cin;
using std::endl;

class MatrixWrongSizeError : public std::logic_error {
public:
    MatrixWrongSizeError() : logic_error("You tried to make operations with different size matrix") {}
};

class MatrixIndexError : public std::out_of_range {
public:
    MatrixIndexError() : out_of_range("You tried to access a nonexistent index") {}
};

// non-specified functions to get "zero" of type T

template<typename T>
T getZero() {
    return T(0);
}

//=============== Matrix class ===============//

template<class T>
class Matrix {
protected:
    int height;
    int width;
    T **matrixTable;

    void setSize(int height, int width);

    void free();

public:
    Matrix(const int height, const int width);

    Matrix(const Matrix<T> &that);

    Matrix<T> &operator=(const Matrix<T> &that);

    int getRowsNumber() const;

    int getColumnsNumber() const;

    Matrix<T> &transpose();

    Matrix<T> getTransposed() const;

    template<class C>
    friend std::ostream &operator<<(std::ostream &os, const Matrix<C> &myMatrix);

    template<class C>
    friend std::istream &operator>>(std::istream &is, Matrix<C> &myMatrix);

    template<class C>
    friend Matrix<C> operator+(const Matrix<C> &first, const Matrix<C> &second);

    template<class C>
    friend Matrix<C> operator-(const Matrix<C> &first, const Matrix<C> &second);

    template<class C>
    friend Matrix<C> operator*(const Matrix<C> &first, const Matrix<C> &second);

    template<class C>
    friend Matrix<C> operator*(const Matrix<C> &myMatrix, const C num);

    template<class C>
    friend Matrix<C> operator*(const C num, const Matrix<C> &myMatrix);

    Matrix<T> &operator+=(const Matrix<T> &that);

    Matrix<T> &operator-=(const Matrix<T> &that);

    Matrix<T> &operator*=(const Matrix<T> &that);

    Matrix<T> &operator*=(const T num);

    T operator()(const int indexI, const int indexJ) const;

    T &operator()(const int indexI, const int indexJ);

    ~Matrix();
};

//================ class Rational ===============//

class RationalDivisionByZero : public std::logic_error {
public:
    RationalDivisionByZero() : logic_error("Division by zero") {}
};


class Rational {
private:
    int p;
    int q;

    int gcd(int x, int y);

    void reduce();

public:
    Rational() : p(0), q(1) {}

    Rational(int p, int q) : p(p), q(q) {
        this->reduce();
    }

    Rational(int p) : p(p), q(1) {}

    friend std::istream &operator>>(std::istream &is, Rational &myRational);

    friend std::ostream &operator<<(std::ostream &os, const Rational &myRational);

    friend Rational operator+(const Rational &first, const Rational &second);

    friend Rational operator-(const Rational &first, const Rational &second);

    friend Rational operator*(const Rational &first, const Rational &second);

    friend Rational operator/(const Rational &first, const Rational &second);

    Rational &operator+=(const Rational &other);

    Rational &operator-=(const Rational &other);

    Rational &operator*=(const Rational &other);

    Rational &operator/=(const Rational &other);

    Rational &operator++();

    Rational operator++(int);

    Rational &operator--();

    Rational operator--(int);

    Rational operator-() const;

    Rational operator+() const;

    friend bool operator<(const Rational &first, const Rational &second);

    friend bool operator>(const Rational &first, const Rational &second);

    friend bool operator<=(const Rational &first, const Rational &second);

    friend bool operator>=(const Rational &first, const Rational &second);

    friend bool operator==(const Rational &first, const Rational &second);

    friend bool operator!=(const Rational &first, const Rational &second);

};


int main() {
    int m, n, p, q;
    cin >> m >> n >> p >> q;

    Matrix<int> A(m, n), B(p, q);
    cin >> A >> B;

    A = A;
    try {
        cout << A + B * 2 - m * A << endl;
        cout << (A -= B += A *= 2) << endl;
        cout << (((A -= B) += A) *= 2) << endl;
    } catch (const MatrixWrongSizeError &) {
        cout << "A and B are of different size." << endl;
    }
    B = A;

    {
        Matrix<int> AA(A);
        Matrix<int> AAA(1, 1);
        AAA = A;
        cout << AA << endl;
        cout << (AAA += Matrix<int>(m, n)) + B << endl;
    }

    Rational r;
    cin >> r;
    Matrix<Rational> C(m, n), D(p, q);
    cin >> C >> D;
    try {
        cout << C * D << endl;
        cout << (C *= D) << endl;
        cout << C << endl;
    } catch (const MatrixWrongSizeError &) {
        cout << "C and D have not appropriate sizes for multiplication." << endl;
    }
    cout << C.getTransposed() * (r * C) << endl;
    cout << C.transpose() << endl;
    try {
        (C(0, 0) *= 6) /= 3;
        cout << C(0, 0) << endl;
        cout << C(m, m) << endl;
    } catch (const MatrixIndexError &) {
        cout << "Index out of range." << endl;
    }

    {
        const Matrix<Rational> &rC = C;
        cout << rC << endl;
        cout << rC.getRowsNumber() << ' ' << rC.getColumnsNumber() << ' ' << rC(0, 0) << endl;
        cout << (C = C) * (Rational(1, 2) * rC).getTransposed() << endl;
    }
    return 0;
}

template<class T>
void Matrix<T>::free() {
    for (int i = 0; i < height; ++i) {
        delete[] matrixTable[i];
    }
    delete[] matrixTable;
    height = 0;
    width = 0;
}

template<class T>
void Matrix<T>::setSize(int height, int width) {
    matrixTable = new T *[height];
    for (int i = 0; i < height; ++i) {
        matrixTable[i] = new T[width];
    }
}

template<class T>
Matrix<T>::Matrix(const int height, const int width) {
    this->height = height;
    this->width = width;
    this->setSize(height, width);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            matrixTable[i][j] = getZero<T>();
        }
    }
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &that) {
    this->height = that.height;
    this->width = that.width;
    this->setSize(height, width);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            matrixTable[i][j] = that.matrixTable[i][j];
        }
    }
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &that) {
    if (this->matrixTable == that.matrixTable) {
        return *this;
    }
    this->free();
    height = that.height;
    width = that.width;
    this->setSize(height, width);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            matrixTable[i][j] = that.matrixTable[i][j];
        }
    }
    return *this;
}

template<class T>
Matrix<T>::~Matrix() {
    this->free();
}

template<class T>
int Matrix<T>::getRowsNumber() const {
    return height;
}

template<class T>
int Matrix<T>::getColumnsNumber() const {
    return width;
}

template<class T>
Matrix<T> &Matrix<T>::transpose() {
    *this = getTransposed();
    return *this;
}

template<class T>
Matrix<T> Matrix<T>::getTransposed() const {
    Matrix<T> temp(this->width, this->height);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            temp.matrixTable[j][i] = this->matrixTable[i][j];
        }
    }
    return temp;
}

template<class T>
std::ostream &operator<<(std::ostream &os, const Matrix<T> &myMatrix) {
    for (int i = 0; i < myMatrix.height; ++i) {
        for (int j = 0; j < myMatrix.width; ++j) {
            os << myMatrix.matrixTable[i][j] << " ";
        }
        os << endl;
    }
    return os;
}

template<class T>
std::istream &operator>>(std::istream &is, Matrix<T> &myMatrix) {
    for (int i = 0; i < myMatrix.height; ++i) {
        for (int j = 0; j < myMatrix.width; ++j) {
            is >> myMatrix.matrixTable[i][j];
        }
    }
    return is;
}

template<class T>
Matrix<T> operator+(const Matrix<T> &first, const Matrix<T> &second) {
    if (first.height != second.height || first.width != second.width) {
        throw MatrixWrongSizeError();
    }
    Matrix<T> sumMatrix(first.height, first.width);
    for (int i = 0; i < sumMatrix.height; ++i) {
        for (int j = 0; j < sumMatrix.width; ++j) {
            sumMatrix.matrixTable[i][j] = first.matrixTable[i][j] + second.matrixTable[i][j];
        }
    }
    return sumMatrix;
}

template<class T>
Matrix<T> operator-(const Matrix<T> &first, const Matrix<T> &second) {
    if (first.height != second.height || first.width != second.width) {
        throw MatrixWrongSizeError();
    }
    Matrix<T> diffMatrix(first.height, first.width);
    for (int i = 0; i < diffMatrix.height; ++i) {
        for (int j = 0; j < diffMatrix.width; ++j) {
            diffMatrix.matrixTable[i][j] = first.matrixTable[i][j] - second.matrixTable[i][j];
        }
    }
    return diffMatrix;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &first, const Matrix<T> &second) {
    if (first.width != second.height) {
        throw MatrixWrongSizeError();
    }
    Matrix<T> multiMatrix(first.height, second.width);
    for (int i = 0; i < first.height; ++i) {
        for (int j = 0; j < second.width; ++j) {
            for (int k = 0; k < first.width; ++k) {
                multiMatrix.matrixTable[i][j] += first.matrixTable[i][k] * second.matrixTable[k][j];
            }
        }
    }
    return multiMatrix;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &myMatrix, const T num) {
    Matrix<T> multiMatrix(myMatrix.height, myMatrix.width);
    for (int i = 0; i < multiMatrix.height; ++i) {
        for (int j = 0; j < multiMatrix.width; ++j) {
            multiMatrix.matrixTable[i][j] = myMatrix.matrixTable[i][j] * num;
        }
    }
    return multiMatrix;
}

template<class T>
Matrix<T> operator*(const T num, const Matrix<T> &myMatrix) {
    return myMatrix * num;
}

template<class T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &that) {
    return *this = (*this + that);
}

template<class T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &that) {
    return *this = (*this - that);
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &that) {
    return *this = (*this * that);
}

template<class T>
T Matrix<T>::operator()(const int indexI, const int indexJ) const {
    if (indexI < 0 || indexI >= this->height || indexJ < 0 || indexJ >= this->width) {
        throw MatrixIndexError();
    }
    return this->matrixTable[indexI][indexJ];
}

template<class T>
T &Matrix<T>::operator()(const int indexI, const int indexJ) {
    if (indexI < 0 || indexI >= this->height || indexJ < 0 || indexJ >= this->width) {
        throw MatrixIndexError();
    }
    return this->matrixTable[indexI][indexJ];
}

template<class T>
Matrix<T> &Matrix<T>::operator*=(const T num) {
    return *this = (*this * num);
}

bool operator<(const Rational &first, const Rational &second) {
    return first.p * second.q < first.q * second.p;
}

bool operator>=(const Rational &first, const Rational &second) {
    return (second < first || first == second);
}

bool operator<=(const Rational &first, const Rational &second) {
    return (first < second || first == second);
}

bool operator>(const Rational &first, const Rational &second) {
    return second < first;
}

bool operator!=(const Rational &first, const Rational &second) {
    return (first < second) || (second < first);
}

bool operator==(const Rational &first, const Rational &second) {
    return !(first < second) && !(second < first);
}

Rational Rational::operator+() const {
    return *this;
}

Rational Rational::operator-() const {
    return Rational(-1 * this->p, this->q);
}

Rational &Rational::operator--() {
    *this -= 1;
    return *this;
}

Rational Rational::operator--(int) {
    Rational copy(this->p, this->q);
    *this -= 1;
    return copy;
}

Rational Rational::operator++(int) {
    Rational copy(this->p, this->q);
    *this += 1;
    return copy;
}

Rational &Rational::operator++() {
    *this += 1;
    return *this;
}

Rational &Rational::operator/=(const Rational &other) {
    *this = *this / other;
    return *this;
}

Rational &Rational::operator*=(const Rational &other) {
    *this = *this * other;
    return *this;
}

Rational &Rational::operator-=(const Rational &other) {
    *this = *this - other;
    return *this;
}

Rational &Rational::operator+=(const Rational &other) {
    *this = *this + other;
    return *this;
}

Rational operator/(const Rational &first, const Rational &second) {
    if (second == 0) {
        throw RationalDivisionByZero();
    }
    return Rational(first.p * second.q, first.q * second.p);
}

Rational operator*(const Rational &first, const Rational &second) {
    return Rational(first.p * second.p, first.q * second.q);
}

Rational operator-(const Rational &first, const Rational &second) {
    return Rational(first.p * second.q - first.q * second.p, first.q * second.q);
}

Rational operator+(const Rational &first, const Rational &second) {
    return Rational(first.p * second.q + first.q * second.p, first.q * second.q);
}

std::ostream &operator<<(std::ostream &os, const Rational &myRational) {
    if (myRational.p == 0 || myRational.q == 1) {
        os << myRational.p;
    } else {
        os << myRational.p << "/" << myRational.q;
    }
    return os;
}

std::istream &operator>>(std::istream &is, Rational &myRational) {
    int temp = scanf("%d/%d", &myRational.p, &myRational.q);
    if (myRational.p == 0 || temp < 2) {
        myRational.q = 1;
    }
    myRational.reduce();
    return is;
}

void Rational::reduce() {
    int cutNumber = gcd(abs(p), abs(q));
    if (q < 0) {
        q *= (-1);
        p *= (-1);
    }
    p /= cutNumber;
    q /= cutNumber;
}

int Rational::gcd(int x, int y) {
    while (y != 0) {
        int c = x % y;
        x = y;
        y = c;
    }
    return x;
}

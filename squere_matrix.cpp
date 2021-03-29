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

class MatrixIsDegenerateError : public std::logic_error {
public:
    MatrixIsDegenerateError() : logic_error(
            "You tried to find the inverse matrix, but it is degenerate. It's impossible") {}
};

// non-specified functions to get "zero" and "one" of type T

template<typename T>
T getZero() {
    return T(0);
}

template<typename T>
T getOne() {
    return T(1);
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

    Matrix<T> &operator*=(const T &num);

    T operator()(const int indexI, const int indexJ) const;

    T &operator()(const int indexI, const int indexJ);

    ~Matrix();
};

//=============== SquareMatrix class ===============//

template<typename T>
class SquareMatrix : public Matrix<T> {
protected:
    template<class C>
    friend void gauss(SquareMatrix<C> &myMatrix, int &counterOfSwaps, bool &zeroDeterminant);

    template<class C>
    friend void makeInvert(SquareMatrix<C> &myMatrix, SquareMatrix<C> &eMatrix, bool &zeroDeterminant);

    template<class C>
    friend void
    special(SquareMatrix<C> &myMatrix, SquareMatrix<C> &eMatrix, bool &zeroDeterminant, bool &isItGauss, int &counter);

public:
    SquareMatrix(const int size);

    SquareMatrix(const Matrix<T> &that);

    SquareMatrix<T> &operator=(const SquareMatrix<T> &that);

    int getSize() const;

    SquareMatrix<T> getInverse() const;

    SquareMatrix<T> &invert();

    SquareMatrix<T> &transpose();

    SquareMatrix<T> getTransposed() const;

    T getTrace() const;

    T getDeterminant() const;

    SquareMatrix<T> &operator+=(const SquareMatrix<T> &that);

    SquareMatrix<T> &operator*=(const SquareMatrix<T> &that);

    SquareMatrix<T> &operator-=(const SquareMatrix<T> &that);

    template<class C>
    friend SquareMatrix<C> operator+(const SquareMatrix<C> &first, const SquareMatrix<C> &second);

    template<class C>
    friend SquareMatrix<C> operator-(const SquareMatrix<C> &first, const SquareMatrix<C> &second);

    template<class C>
    friend SquareMatrix<C> operator*(const SquareMatrix<C> &first, const SquareMatrix<C> &second);

    template<class C>
    friend SquareMatrix<C> operator*(const SquareMatrix<C> &myMatrix, const C &num);

    template<class C>
    friend SquareMatrix<C> operator*(const C &num, const SquareMatrix<C> &myMatrix);
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

    int getNumerator() const;

    int getDenominator() const;

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
    int m, n, p;
    Rational r;
    cin >> m >> n >> p >> r;

    Matrix<Rational> A(m, n);
    SquareMatrix<Rational> S(p);
    cin >> A >> S;

    try {
        cout << (A * S) * A.getTransposed() << endl;
    } catch (const MatrixWrongSizeError &) {
        cout << "A and S have not appropriate sizes for multiplication." << endl;
    }

    cout << (r * (S = S) * S).getSize() << endl;

    SquareMatrix<Rational> P(S);

    cout << (P * (S + S - Rational(3) * P)).getDeterminant() << endl;

    const SquareMatrix<Rational> &rS = S;

    cout << rS.getSize() << ' ' << rS.getDeterminant() << ' ' << rS.getTrace() << endl;
    cout << (S = S) * (S + rS) << endl;
    cout << (S *= S) << endl;

    try {
        cout << rS.getInverse() << endl;
        cout << P.invert().getTransposed().getDeterminant() << endl;
        cout << P << endl;
    } catch (const MatrixIsDegenerateError &) {
        cout << "Cannot inverse matrix." << endl;
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
Matrix<T> &Matrix<T>::operator*=(const T &num) {
    return *this = (*this * num);
}

template<class T>
SquareMatrix<T>::SquareMatrix(const int size) : Matrix<T>(size, size) {}

template<class T>
SquareMatrix<T>::SquareMatrix(const Matrix<T> &that) : Matrix<T>(that) {}

template<class T>
SquareMatrix<T> &SquareMatrix<T>::operator=(const SquareMatrix<T> &that) {
    Matrix<T>::operator=(that);
    return *this;
}

template<class T>
int SquareMatrix<T>::getSize() const {
    return this->width;
}

template<class T>
SquareMatrix<T> operator+(const SquareMatrix<T> &first, const SquareMatrix<T> &second) {
    SquareMatrix<T> sumMatrix = (static_cast<Matrix<T>>(first) + static_cast<Matrix<T>>(second));
    return sumMatrix;
}

template<class T>
SquareMatrix<T> operator-(const SquareMatrix<T> &first, const SquareMatrix<T> &second) {
    SquareMatrix<T> diffMatrix = (static_cast<Matrix<T>>(first) - static_cast<Matrix<T>>(second));
    return diffMatrix;
}

template<class T>
SquareMatrix<T> operator*(const SquareMatrix<T> &first, const SquareMatrix<T> &second) {
    SquareMatrix<T> multiMatrix = (static_cast<Matrix<T>>(first) * static_cast<Matrix<T>>(second));
    return multiMatrix;
}

template<class T>
SquareMatrix<T> operator*(const SquareMatrix<T> &myMatrix, const T &num) {
    const Matrix<T> &temp = myMatrix;
    return temp * num;
}

template<class T>
SquareMatrix<T> operator*(const T &num, const SquareMatrix<T> &myMatrix) {
    return myMatrix * num;
}

template<class T>
T SquareMatrix<T>::getTrace() const {
    T trace = getZero<T>();
    for (int i = 0; i < this->width; ++i) {
        trace += this->matrixTable[i][i];
    }
    return trace;
}

template<class T>
void special(SquareMatrix<T> &myMatrix, SquareMatrix<T> &eMatrix,
             bool &zeroDeterminant, bool &isItGauss, int &counter) {
    int size = myMatrix.getSize();
    bool zeroSubMatrix;
    for (int i = 0; i < size; ++i) {
        if (myMatrix(i, i) != 0) {
            for (int j = i + 1; j < size; ++j) {
                if (myMatrix(j, i) != 0) {
                    T temp = myMatrix(j, i) / myMatrix(i, i);
                    for (int k = 0; k < size; ++k) {
                        myMatrix(j, k) -= myMatrix(i, k) * temp;
                        if (!isItGauss) {
                            eMatrix(j, k) -= eMatrix(i, k) * temp;
                        }
                    }
                }
            }
        } else {
            zeroSubMatrix = true;
            for (int j = i + 1; j < size; ++j) {
                if (myMatrix(j, i) != 0) {
                    std::swap(myMatrix.matrixTable[i], myMatrix.matrixTable[j]);
                    if (!isItGauss) {
                        std::swap(eMatrix.matrixTable[i], eMatrix.matrixTable[j]);
                    } else {
                        ++counter;
                    }
                    zeroSubMatrix = false;
                    break;
                }
            }
            if (i < size - 1) {
                --i;
            }
            if (zeroSubMatrix) {
                zeroDeterminant = true;
                break;
            }
        }
    }
}

template<class T>
void gauss(SquareMatrix<T> &myMatrix, int &counterOfSwaps, bool &zeroDeterminant) {
    bool isItGauss = true;
    SquareMatrix<T> eMatrix(myMatrix.getSize());
    for (int i = 0; i < eMatrix.getSize(); ++i) {
        eMatrix(i, i) = getOne<T>();
    }
    special(myMatrix, eMatrix, zeroDeterminant, isItGauss, counterOfSwaps);
}

template<class T>
T SquareMatrix<T>::getDeterminant() const {
    int size = this->getSize();
    T determinant = getOne<T>();
    int counter = 0;
    SquareMatrix<T> temp(*this);
    bool zeroDeterminant = false;
    gauss(temp, counter, zeroDeterminant);
    if (zeroDeterminant) {
        return getZero<T>();
    } else {
        for (int i = 0; i < size; ++i) {
            determinant *= temp(i, i);
        }
        return determinant * ((counter % 2 == 0) ? 1 : (-1));
    }
}

template<class T>
void makeInvert(SquareMatrix<T> &myMatrix, SquareMatrix<T> &eMatrix, bool &zeroDeterminant) {
    int size = myMatrix.getSize();
    bool isItGauss = false;
    int counter = 0;
    special(myMatrix, eMatrix, zeroDeterminant, isItGauss, counter);
    if (!zeroDeterminant) {
        for (int i = size - 1; i >= 0; --i) {
            for (int j = i - 1; j >= 0; --j) {
                T temp = myMatrix(j, i) / myMatrix(i, i);
                for (int k = size - 1; k >= 0; --k) {
                    myMatrix(j, k) -= myMatrix(i, k) * temp;
                    eMatrix(j, k) -= eMatrix(i, k) * temp;
                }
            }
        }

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                eMatrix(i, j) /= myMatrix(i, i);
            }
            myMatrix(i, i) = 1;
        }
    }
}

template<typename T>
SquareMatrix<T> SquareMatrix<T>::getInverse() const {
    bool zeroDeterminant = false;
    SquareMatrix<T> myMatrix(*this);
    SquareMatrix<T> eMatrix(this->getSize());
    for (int i = 0; i < eMatrix.getSize(); ++i) {
        eMatrix(i, i) = getOne<T>();
    }
    makeInvert(myMatrix, eMatrix, zeroDeterminant);
    if (zeroDeterminant) {
        throw MatrixIsDegenerateError();
    } else {
        return eMatrix;
    }
}


template<typename T>
SquareMatrix<T> &SquareMatrix<T>::invert() {
    *this = getInverse();
    return *this;
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::transpose() {
    *this = this->getTransposed();
    return *this;
}

template<typename T>
SquareMatrix<T> SquareMatrix<T>::getTransposed() const {
    return Matrix<T>::getTransposed();
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::operator+=(const SquareMatrix<T> &that) {
    return *this = (*this + that);
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::operator*=(const SquareMatrix<T> &that) {
    return *this = (*this * that);
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::operator-=(const SquareMatrix<T> &that) {
    return *this = (*this - that);
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
    if (myRational.p == 0 or myRational.q == 1) {
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

int Rational::getDenominator() const {
    return q;
}

int Rational::getNumerator() const {
    return p;
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
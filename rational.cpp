#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <cstdlib>

using std::cout;
using std::cin;
using std::endl;

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
    int a;
    cin >> a;

    int p, q;
    cin >> p >> q;
    const Rational rc(p, q);
    cout << rc.getNumerator() << ' ' << rc.getDenominator() << endl;

    Rational r1, r2;
    cin >> r1 >> r2;

    cout << r1 << endl;
    cout << r2 << endl;

    try {
        cout << 1 / r1 << endl;
    } catch (const RationalDivisionByZero &ex) {
        cout << "Cannot get reciprocal of r1." << endl;
    }

    try {
        cout << rc / r2 << endl;
    } catch (const RationalDivisionByZero &ex) {
        cout << "Cannot divide by r2." << endl;
    }

    cout << (r1 < r2) << endl;
    cout << (r1 <= r2) << endl;
    cout << (r1 > r2) << endl;
    cout << (r1 >= r2) << endl;
    cout << (r1 == r2) << endl;
    cout << (r1 != r2) << endl;

    cout << (r1 < a) << endl;
    cout << (r1 <= a) << endl;
    cout << (r1 > a) << endl;
    cout << (r1 >= a) << endl;
    cout << (r1 == a) << endl;
    cout << (r1 != a) << endl;

    cout << (a < r2) << endl;
    cout << (a <= r2) << endl;
    cout << (a > r2) << endl;
    cout << (a >= r2) << endl;
    cout << (a == r2) << endl;
    cout << (a != r2) << endl;

    cout << rc + a << endl
         << a + rc << endl
         << -rc * r1 << endl
         << (+r1 - r2 * rc) * a << endl;

    cout << ++r1 << endl;
    cout << r1 << endl;
    cout << r1++ << endl;
    cout << r1 << endl;
    cout << --r1 << endl;
    cout << r1 << endl;
    cout << r1-- << endl;
    cout << r1 << endl;
    cout << ++ ++r1 << endl;
    cout << r1 << endl;

    cout << ((((r1 += r2) /= Rational(-5, 3)) -= rc) *= a) << endl;
    cout << (r1 += r2 /= 3) << endl;
    cout << r1 << endl;
    cout << r2 << endl;
    return 0;
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

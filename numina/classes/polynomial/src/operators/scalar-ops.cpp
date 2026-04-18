#include "../../include/numina/polynomial.h"
// Created by Vadim on 18.06.2025.
Polynomial& Polynomial::operator+=(const Type& scalar) {
    c[n - 1] += scalar;
    return *this;
}

Polynomial& Polynomial::operator-=(const Type& scalar) {
    c[n - 1] -= scalar;
    return *this;
}

Polynomial& Polynomial::operator*=(const Type& scalar) {
    if (scalar == Type(0)) {
        zeroing();
        return *this;
    }

    for (std::size_t i = 0; i < n; ++i)
        c[i] *= scalar;
    return *this;
}

Polynomial& Polynomial::operator/=(const Type& scalar) {
    for (std::size_t i = 0; i < n; ++i)
        c[i] /= scalar;
    return *this;
}

Polynomial operator+(const Polynomial& poly, const Polynomial::Type& scalar) {
    auto p = poly;
    return p += scalar;
}

Polynomial operator-(const Polynomial& poly, const Polynomial::Type& scalar) {
    auto p = poly;
    return p -= scalar;
}

Polynomial operator*(const Polynomial& poly, const Polynomial::Type& scalar) {
    auto p = poly;
    return p *= scalar;
}

Polynomial operator/(const Polynomial& poly, const Polynomial::Type& scalar) {
    auto p = poly;
    return p /= scalar;
}

Polynomial operator+(const Polynomial::Type& scalar, const Polynomial& poly) {
    auto p = poly;
    return p += scalar;
}

Polynomial operator-(const Polynomial::Type& scalar, const Polynomial& poly) {
    auto p = -poly;
    return p += scalar;
}

Polynomial operator*(const Polynomial::Type& scalar, const Polynomial& poly) {
    auto p = poly;
    return p *= scalar;
}
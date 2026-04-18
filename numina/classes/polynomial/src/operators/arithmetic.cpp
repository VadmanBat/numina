#include "numina/polynomial.h"
// Created by Vadim on 18.06.2025.
Polynomial Polynomial::operator+() const {
    return {coeffs};
}

Polynomial Polynomial::operator-() const {
    std::vector<Type> answer(n);
    const auto a = answer.data();
    for (std::size_t i = 0; i < n; ++i)
        a[i] = -c[i];
    return {std::move(answer)};
}

Polynomial Polynomial::operator+(const Polynomial& other) const {
    std::vector<Type> answer;

    if (n < other.n) {
        answer       = other.coeffs;
        const auto a = answer.data() + other.n - n;
        for (std::size_t i = 0; i < n; ++i)
            a[i] += c[i];
        return {std::move(answer)};
    }

    answer       = coeffs;
    const auto a = answer.data() + n - other.n;
    for (std::size_t i = 0; i < other.n; ++i)
        a[i] += other.c[i];
    return {std::move(answer)};
}

Polynomial Polynomial::operator-(const Polynomial& other) const {
    std::vector<Type> answer;

    if (n < other.n) {
        answer.assign(other.c, other.c + other.n);
        auto a       = answer.data();
        const auto d = other.n - n;
        for (std::size_t i = 0; i < d; ++i)
            a[i] = -a[i];
        a += d;
        for (std::size_t i = 0; i < n; ++i)
            a[i] = c[i] - a[i];
        return {std::move(answer)};
    }

    answer       = coeffs;
    const auto a = answer.data() + n - other.n;
    for (std::size_t i = 0; i < other.n; ++i)
        a[i] -= other.c[i];
    return {std::move(answer)};
}

Polynomial Polynomial::operator*(const Polynomial& other) const {
    if (isZero() || other.isZero())
        return {};

    std::vector<Type> answer(n + other.n - 1, 0);
    const auto a = answer.data();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < other.n; ++j)
            a[i + j] += c[i] * other.c[j];
    return {std::move(answer)};
}

Polynomial& Polynomial::operator+=(const Polynomial& other) {
    if (n < other.n) {
        std::vector<Type> answer = other.coeffs;
        const auto a             = answer.data() + other.n - n;
        for (std::size_t i = 0; i < n; ++i)
            a[i] += c[i];
        return *this = std::move(answer);
    }

    const auto a = c + n - other.n;
    for (std::size_t i = 0; i < other.n; ++i)
        a[i] += other.c[i];

    remove_leading_zeros();
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& other) {
    if (n < other.n) {
        std::vector<Type> answer = other.coeffs;
        auto a                   = answer.data();
        const auto d             = other.n - n;
        for (std::size_t i = 0; i < d; ++i)
            a[i] = -a[i];
        a += d;
        for (std::size_t i = 0; i < n; ++i)
            a[i] = c[i] - a[i];
        return *this = std::move(answer);
    }

    const auto a = c + n - other.n;
    for (std::size_t i = 0; i < other.n; ++i)
        a[i] -= other.c[i];

    remove_leading_zeros();
    return *this;
}

Polynomial& Polynomial::operator*=(const Polynomial& other) {
    if (isZero())
        return *this;
    if (other.isZero()) {
        zeroing();
        return *this;
    }

    std::vector<Type> answer(n + other.n - 1, 0);
    const auto a = answer.data();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < other.n; ++j)
            a[i + j] += c[i] * other.c[j];
    return *this = std::move(answer);
}

Polynomial& Polynomial::operator/=(const Polynomial& divisor) {
    return *this = *this / divisor;
}

Polynomial& Polynomial::operator%=(const Polynomial& divisor) {
    return *this = *this % divisor;
}
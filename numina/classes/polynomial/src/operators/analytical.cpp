#include "numina/polynomial.h"
// Created by Vadim on 18.06.2025.
Polynomial::Type Polynomial::operator()(const Type& x) const {
    Type result = c[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c[i];
    return result;
}

Polynomial::Complex Polynomial::operator()(const Complex& x) const {
    Complex result = c[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c[i];
    return result;
}

Polynomial::LongType Polynomial::operator()(const LongType& x) const {
    LongType result = c[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c[i];
    return result;
}

Polynomial::LongComplex Polynomial::operator()(const LongComplex& x) const {
    LongComplex result = c[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c[i];
    return result;
}

Polynomial Polynomial::derivative() const {
    if (n < 2) /// n = 1: C -> 0
        return {};

    std::vector<Type> new_coeffs(n - 1);
    const auto new_c = new_coeffs.data();
    int power        = static_cast<int>(n);
    for (std::size_t i = 0; i < n - 1; ++i)
        new_c[i] = --power * c[i];

    return {std::move(new_coeffs)};
}

Polynomial Polynomial::derivative(const std::size_t k) const {
    if (k == 0)
        return *this;
    if (n <= k)
        return {};

    const auto new_n = n - k;
    std::vector<Type> new_coeffs(new_n);
    const auto new_c = new_coeffs.data();
    auto power       = n - 1;
    auto mul         = power;
    for (std::size_t j = 1; j < k; ++j)
        mul *= --power;
    for (std::size_t i = 0; i < new_n; ++i) {
        new_c[i] = static_cast<Type>(mul) * c[i];
        (mul     /= --power + k) *= power;
    }

    return {std::move(new_coeffs)};
}

void Polynomial::makeDerivative() {
    if (n < 2) {
        zeroing();
        return;
    }

    int power = static_cast<int>(n);
    --n;
    for (std::size_t i = 0; i < n; ++i)
        c[i] *= --power;
    coeffs.pop_back();
}

void Polynomial::makeDerivative(const std::size_t k) {
    if (k == 0)
        return;
    if (n <= k) {
        zeroing();
        return;
    }

    auto power = n - 1;
    n          -= k;
    auto mul   = power;
    for (std::size_t j = 1; j < k; ++j)
        mul *= --power;
    for (std::size_t i = 0; i < n; ++i) {
        c[i] *= static_cast<Type>(mul);
        (mul /= --power + k) *= power;
    }
    coeffs.resize(n);
}

Polynomial Polynomial::integral(const Type& constant) const {
    std::vector<Type> new_coeffs(n + 1);
    const auto new_c = new_coeffs.data();
    int power        = static_cast<int>(n);
    for (std::size_t i = 0; i < n; ++i) {
        new_c[i] = c[i] / power;
        --power;
    }
    new_c[n] = constant;
    return {std::move(new_coeffs)};
}

void Polynomial::makeIntegral(const Type& constant) {
    int power = static_cast<int>(n);
    for (std::size_t i = 0; i < n; ++i) {
        c[i] /= power;
        --power;
    }
    ++n;
    coeffs.push_back(constant);
    c = coeffs.data();
}

Polynomial Polynomial::deflate(const Type& root) const {
    if (n < 3) /// n = 2: x + 1 -> 0
        return {};

    const auto new_n = n - 1;
    std::vector<Type> answer(new_n);
    const auto a = answer.data();

    a[0] = c[0];
    for (std::size_t i = 1; i < new_n; ++i)
        a[i] = c[i] + root * a[i - 1];

    return {std::move(answer)};
}

void Polynomial::makeDeflate(const Type& root) {
    if (n < 3) {
        zeroing();
        return;
    }

    --n;
    for (std::size_t i = 1; i < n; ++i)
        c[i] += root * c[i - 1];
    coeffs.pop_back();
}

Polynomial Polynomial::deflateConjRoot(const Complex& root) const {
    if (n < 4) /// n = 3: x^2 + x + 1 -> 0
        return {};

    const auto new_n = n - 2;

    const Type real = root.real();
    const Type imag = root.imag();
    const Type c1   = 2 * real;
    const Type c2   = real * real + imag * imag;

    std::vector<Type> answer(new_n);
    const auto a = answer.data();

    a[0] = c[0];
    a[1] = c[1] + c1 * a[0];
    for (std::size_t i = 2; i < new_n; ++i)
        a[i] = c[i] + c1 * a[i - 1] - c2 * a[i - 2];

    return {std::move(answer)};
}

void Polynomial::makeDeflateConjRoot(const Complex& root) {
    if (n < 4) {
        zeroing();
        return;
    }

    const Type real = root.real();
    const Type imag = root.imag();
    const Type c1   = 2 * real;
    const Type c2   = real * real + imag * imag;

    n    -= 2;
    c[1] += c1 * c[0];
    for (std::size_t i = 2; i < n; ++i)
        c[i] += c1 * c[i - 1] - c2 * c[i - 2];
    coeffs.resize(n);
}

Polynomial Polynomial::monic() const {
    if (c[0] == 1)
        return {coeffs};

    std::vector<Type> answer(coeffs);
    const auto a = answer.data();
    for (std::size_t i = 1; i < n; ++i)
        a[i] /= a[0];
    a[0] = 1;

    return {std::move(answer)};
}

void Polynomial::makeMonic() {
    if (c[0] == 1)
        return;

    for (std::size_t i = 1; i < n; ++i)
        c[i] /= c[0];
    c[0] = 1;
}
#include "numina/polynomial.h"
// Created by Vadim on 21.08.2025.
Polynomial::Type Polynomial::multiplicity(const Type& x, const std::size_t m) const {
    const auto ax = derivative(m - 1)(x);
    auto bx2      = derivative(m)(x);
    bx2           *= bx2;
    const auto cx = derivative(m + 1)(x);
    return static_cast<Type>(m) - 1 + bx2 / (bx2 - ax * cx);
}

Polynomial::Type Polynomial::multiplicity(const Complex& x, const std::size_t m) const {
    const auto ax = derivative(m - 1)(x);
    auto bx2      = derivative(m)(x);
    bx2           *= bx2;
    const auto cx = derivative(m + 1)(x);
    return static_cast<Type>(m) - 1 + (bx2 / (bx2 - ax * cx)).real();
}

std::size_t Polynomial::computeMultiplicity(const Type& x) const {
    std::vector<Type> dfx_1;
    std::vector<LongType> dfx_2;
    dfx_1.reserve(n);
    dfx_2.reserve(n);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        dfx_1.push_back(deriv(static_cast<Type>(x)));
        dfx_2.push_back(deriv(static_cast<LongType>(x)));
        deriv.makeDerivative();
    }

    int prev_m = -1;
    for (std::size_t m = 1; m <= n; ++m) {
        const Type m1 = static_cast<Type>(m) - 1 + dfx_1[m] * dfx_1[m] / (
                            dfx_1[m] * dfx_1[m] - dfx_1[m - 1] * dfx_1[m + 1]);
        const LongType m2 = static_cast<LongType>(m) - 1 + dfx_2[m] * dfx_2[m] / (
                                dfx_2[m] * dfx_2[m] - dfx_2[m - 1] * dfx_2[m + 1]);

        if (std::abs(m1 - m2) < 1e-6)
            if (m == std::round(m2))
                return m;
        if (const auto rounded = static_cast<int>(std::round(m1)); rounded == std::round(m2)) {
            if (prev_m >= static_cast<int>(m) && prev_m == rounded)
                return prev_m;
            prev_m = rounded;
        }

        dfx_1.push_back(deriv(static_cast<Type>(x)));
        dfx_2.push_back(deriv(static_cast<LongType>(x)));
        deriv.makeDerivative();
    }

    return 1;
}

std::size_t Polynomial::computeMultiplicity(const Complex& x) const {
    std::vector<Complex> dfx_1;
    std::vector<LongComplex> dfx_2;
    dfx_1.reserve(n);
    dfx_2.reserve(n);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        dfx_1.push_back(deriv(static_cast<Complex>(x)));
        dfx_2.push_back(deriv(static_cast<LongComplex>(x)));
        deriv.makeDerivative();
    }

    int prev_m = -1;
    for (std::size_t m = 1; m <= n; ++m) {
        const Type m1 = static_cast<Type>(m) - 1 + (dfx_1[m] * dfx_1[m] / (
                               dfx_1[m] * dfx_1[m] - dfx_1[m - 1] * dfx_1[m + 1])).real();
        const LongType m2 = static_cast<LongType>(m) - 1 + (dfx_2[m] * dfx_2[m] / (
                                   dfx_2[m] * dfx_2[m] - dfx_2[m - 1] * dfx_2[m + 1])).real();

        if (std::abs(m1 - m2) < 1e-6)
            if (m == std::round(m2))
                return m;
        if (const auto rounded = static_cast<int>(std::round(m1)); rounded == std::round(m2)) {
            if (prev_m >= static_cast<int>(m) && prev_m == rounded)
                return prev_m;
            prev_m = rounded;
        }

        dfx_1.push_back(deriv(static_cast<Complex>(x)));
        dfx_2.push_back(deriv(static_cast<LongComplex>(x)));
        deriv.makeDerivative();
    }

    return 1;
}

/*
#include "numina/quadratic-deriv-poly.h"

template <typename T>
struct DerivativeValues {
    std::vector<T> df; // все производные в точке x
    T num = 0;         // значение числителя
    T den = 0;         // значение знаменателя
    T m   = 0;         // оценка кратности на текущей итерации

    DerivativeValues() = default;

    explicit DerivativeValues(std::size_t reserve_size) {
        df.reserve(reserve_size);
    }
};

std::size_t Polynomial::computeMultiplicity(const Type& x) const {
    numina::QuadraticDerivPoly num{{1, 1, 1}}, den{{1, 1, 1}, {0, 2, -1}};

    DerivativeValues<Type> t1(n);
    DerivativeValues<LongType> t2(n);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        t1.df.push_back(deriv(static_cast<Type>(x)));
        t2.df.push_back(deriv(static_cast<LongType>(x)));
        deriv.makeDerivative();
    }

    int prev_m = -1;
    for (std::size_t m = 1; m <= n; ++m) {
        t1.num = t1.den = 0;
        t2.num = t2.den = 0;

        for (const auto [pr, coeff] : num.getTerms()) {
            const auto [a, b] = pr;
            t1.num            += static_cast<Type>(coeff) * t1.df[a] * t1.df[b];
            t2.num            += static_cast<LongType>(coeff) * t2.df[a] * t2.df[b];
        }
        for (const auto [pr, coeff] : den.getTerms()) {
            const auto [a, b] = pr;
            t1.den            += static_cast<Type>(coeff) * t1.df[a] * t1.df[b];
            t2.den            += static_cast<LongType>(coeff) * t2.df[a] * t2.df[b];
        }

        t1.m = t1.num / t1.den;
        t2.m = t2.num / t2.den;

        if (std::abs(t1.m - t2.m) < 1e-6)
            if (m == std::round(t2.m))
                return m;
        if (const auto rounded = static_cast<int>(std::round(t1.m)); rounded == std::round(t2.m)) {
            if (prev_m >= static_cast<int>(m) && prev_m == rounded)
                return prev_m;
            prev_m = rounded;
        }

        num.doubleDifferentiate();
        den.doubleDifferentiate();
        for (int i = 0; i < 2; ++i) {
            t1.df.push_back(deriv(static_cast<Type>(x)));
            t2.df.push_back(deriv(static_cast<LongType>(x)));
            deriv.makeDerivative();
        }
    }

    return 1;
}*/
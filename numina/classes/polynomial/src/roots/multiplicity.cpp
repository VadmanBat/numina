//
// Created by Vadim on 21.08.2025.
//

#include "../../include/numina/polynomial.h"
#include "numina/quadratic-deriv-poly.h"
#include <iostream>

Polynomial::Type Polynomial::multiplicity(const Type& x) const {
    const auto d1f = derivative();
    const auto d2f = d1f.derivative();
    auto s         = d1f(x);
    s              *= s;
    return std::abs(s / (s - (*this)(x) * d2f(x)));
}

Polynomial::Type Polynomial::multiplicity(const Complex& x) const {
    const auto d1f = derivative();
    const auto d2f = d1f.derivative();
    auto s         = d1f(x);
    s              *= s;
    return std::abs(s / (s - (*this)(x) * d2f(x)));
}

Polynomial::Type Polynomial::multiplicity(const Type& x, int m) const {
    // Предвычисляем производные (достаточно для полинома 10-й степени)
    const auto d1f  = derivative();     // p'
    const auto d2f  = d1f.derivative(); // p''
    const auto d3f  = d2f.derivative(); // p'''
    const auto d4f  = d3f.derivative(); // p⁽⁴⁾
    const auto d5f  = d4f.derivative(); // p⁽⁵⁾
    const auto d6f  = d5f.derivative(); // p⁽⁶⁾
    const auto d7f  = d6f.derivative(); // p⁽⁷⁾
    const auto d8f  = d7f.derivative(); // p⁽⁸⁾
    const auto d9f  = d8f.derivative(); // p⁽⁹⁾
    const auto d10f = d9f.derivative(); // p⁽¹⁰⁾

    //std::cout << (*this)(x) << ", " << d2f(x) << ", " << d1f(x) << '\n';

    switch (m) {
        case 1: // 0 применений Лопиталя (твоя исходная переписанная форма)
            return Type(1) / (Type(1) - (*this)(x) * d2f(x) / (d1f(x) * d1f(x)));

        case 2: // После 2 применений Лопиталя
            return Type(2) * (d2f(x) * d2f(x) + d1f(x) * d3f(x))
                   / (d2f(x) * d2f(x) - (*this)(x) * d4f(x));

        case 3: // После 3 применений Лопиталя
            return Type(2) * (d1f(x) * d5f(x) + Type(4) * d2f(x) * d4f(x) + Type(3) * d3f(x) * d3f(x))
                   / (-(*this)(x) * d6f(x) - Type(2) * d1f(x) * d5f(x) + d2f(x) * d4f(x) + Type(2) * d3f(x) * d3f(x));

        case 4: // После 5 применений Лопиталя
            return Type(2) * (-d1f(x) * d7f(x) - Type(6) * d2f(x) * d6f(x) - Type(15) * d3f(x) * d5f(x) - Type(10) *
                              d4f(x) * d4f(x))
                   / ((*this)(x) * d8f(x) + Type(4) * d1f(x) * d7f(x) + Type(4) * d2f(x) * d6f(x)
                      - Type(4) * d3f(x) * d5f(x) - Type(5) * d4f(x) * d4f(x));

        case 5:
            return Type(2) * (-d1f(x) * d9f(x) - Type(8) * d2f(x) * d8f(x)
                              - Type(28) * d3f(x) * d7f(x) - Type(56) * d4f(x) * d6f(x)
                              - Type(35) * d5f(x) * d5f(x))
                   / ((*this)(x) * d10f(x) + Type(6) * d1f(x) * d9f(x) + Type(13) * d2f(x) * d8f(x)
                      + Type(8) * d3f(x) * d7f(x) - Type(14) * d4f(x) * d6f(x)
                      - Type(14) * d5f(x) * d5f(x));

        default:
            return Type(0); // или throw std::invalid_argument("m > 6 не поддерживается");
    }
}

#include <iostream>

int Polynomial::__computeMultiplicity(const Type& x, bool d) const {
    numina::QuadraticDerivPoly num{{1, 1, 1}}, den{{1, 1, 1}, {0, 2, -1}};
    std::vector<Type> df_1;
    std::vector<long double> df_2;
    df_1.reserve(n);
    df_2.reserve(n);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        df_1.push_back(deriv(static_cast<Type>(x)));
        df_2.push_back(deriv(static_cast<long double>(x)));
        deriv.makeDerivative();
    }

    int prev_m = -1;
    for (int m = 1; m <= n; ++m) {
        Type num_val_1        = 0, den_val_1 = 0;
        long double num_val_2 = 0, den_val_2 = 0;

        for (const auto [pr, c] : num.getTerms()) {
            const auto [a, b] = pr;
            num_val_1         += static_cast<Type>(c) * df_1[a] * df_1[b];
            num_val_2         += static_cast<long double>(c) * df_2[a] * df_2[b];
        }
        for (const auto [pr, c] : den.getTerms()) {
            const auto [a, b] = pr;
            den_val_1         += static_cast<Type>(c) * df_1[a] * df_1[b];
            den_val_2         += static_cast<long double>(c) * df_2[a] * df_2[b];
        }

        const auto m1 = num_val_1 / den_val_1;
        const auto m2 = num_val_2 / den_val_2;
        if (d)
            std::cout << "f1: " << (std::abs(m1 - m2) < 1e-6) << ", i = " << m << '\n';
        if (std::abs(m1 - m2) < 1e-6)
            if (m == std::round(m2))
                return m;

        if (std::round(m1) == std::round(m2)) {
            if (prev_m >= m && prev_m == std::round(m1)) {
                if (d)
                    std::cout << "f2:\n";
                return prev_m;
            }
            prev_m = std::round(m1);
        }

        num.doubleDifferentiate();
        den.doubleDifferentiate();

        for (int i = 0; i < 2; ++i) {
            df_1.push_back(deriv(static_cast<Type>(x)));
            df_2.push_back(deriv(static_cast<long double>(x)));
            deriv.makeDerivative();
        }
    }
    return 1;
}

int Polynomial::__computeMultiplicity(const Complex& x, bool d) const {
    numina::QuadraticDerivPoly num{{1, 1, 1}}, den{{1, 1, 1}, {0, 2, -1}};
    std::vector<Complex> df_1;
    std::vector<LongComplex> df_2;
    df_1.reserve(n);
    df_2.reserve(n);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        df_1.push_back(deriv(static_cast<Complex>(x)));
        df_2.push_back(deriv(static_cast<LongComplex>(x)));
        deriv.makeDerivative();
    }

    int prev_m = -1;
    for (int m = 1; m <= n; ++m) {
        Complex num_val_1     = 0, den_val_1 = 0;
        LongComplex num_val_2 = 0, den_val_2 = 0;

        for (const auto [pr, c] : num.getTerms()) {
            const auto [a, b] = pr;
            num_val_1         += static_cast<Type>(c) * df_1[a] * df_1[b];
            num_val_2         += static_cast<long double>(c) * df_2[a] * df_2[b];
        }
        for (const auto [pr, c] : den.getTerms()) {
            const auto [a, b] = pr;
            den_val_1         += static_cast<Type>(c) * df_1[a] * df_1[b];
            den_val_2         += static_cast<long double>(c) * df_2[a] * df_2[b];
        }

        const auto m1 = std::abs(num_val_1 / den_val_1);
        const auto m2 = std::abs(num_val_2 / den_val_2);
        //std::cout << m1 << ' ' << m2 << '\n';
        if (d)
            std::cout << "f1: " << (std::abs(m1 - m2) < 1e-6) << ", i = " << m << '\n';
        if (std::abs(m1 - m2) < 1e-6)
            if (m == std::round(m2))
                return m;

        if (std::round(m1) == std::round(m2)) {
            if (prev_m >= m && prev_m == std::round(m2)) {
                if (d)
                    std::cout << "f2:\n";
                return prev_m;
            }
            prev_m = static_cast<int>(std::round(m2));
        }

        num.doubleDifferentiate();
        den.doubleDifferentiate();

        for (int i = 0; i < 2; ++i) {
            df_1.push_back(deriv(static_cast<Complex>(x)));
            df_2.push_back(deriv(static_cast<LongComplex>(x)));
            deriv.makeDerivative();
        }
    }
    return 1;
}

double Polynomial::lag(const Type& x, int m) const {
    const auto a = derivative(m - 1);
    const auto b = derivative(m);
    const auto c = derivative(m + 1);
    return m - 1 + b(x) * b(x) / (b(x) * b(x) - a(x) * c(x));
}

long double Polynomial::lag(const long double& x, int m) const {
    const auto a = derivative(m - 1);
    const auto b = derivative(m);
    const auto c = derivative(m + 1);
    return m - 1 + b(x) * b(x) / (b(x) * b(x) - a(x) * c(x));
}

long double Polynomial::multiplicity(const long double& x, int m) const {
    // Предвычисляем производные (достаточно для полинома 10-й степени)

    auto d0f  = *this;            // p
    auto d1f  = d0f.derivative(); // p'
    auto d2f  = d1f.derivative(); // p''
    auto d3f  = d2f.derivative(); // p'''
    auto d4f  = d3f.derivative(); // p⁽⁴⁾
    auto d5f  = d4f.derivative(); // p⁽⁵⁾
    auto d6f  = d5f.derivative(); // p⁽⁶⁾
    auto d7f  = d6f.derivative(); // p⁽⁷⁾
    auto d8f  = d7f.derivative(); // p⁽⁸⁾
    auto d9f  = d8f.derivative(); // p⁽⁹⁾
    auto d10f = d9f.derivative(); // p⁽¹⁰⁾

    long double v[11] = {
        (*this)(x),
        d1f(x), d2f(x), d3f(x), d4f(x), d5f(x), d6f(x), d7f(x), d8f(x), d9f(x), d10f(x),
    };

    double y      = x;
    double dv[11] = {
        (*this)(y),
        d1f(y), d2f(y), d3f(y), d4f(y), d5f(y), d6f(y), d7f(y), d8f(y), d9f(y), d10f(y),
    };

    switch (m) {
        case 1: // 0 применений Лопиталя (твоя исходная переписанная форма)
            return (long double)(1) / ((long double)(1) - d0f(x) * d2f(x) / (d1f(x) * d1f(x)));

        case 2: // После 2 применений Лопиталя
            return (long double)(2) * (d2f(x) * d2f(x) + d1f(x) * d3f(x))
                   / (d2f(x) * d2f(x) - d0f(x) * d4f(x));

        case 3: // После 3 применений Лопиталя
            return (long double)(2) * (d1f(x) * d5f(x) + (long double)(4) * d2f(x) * d4f(x) + (long double)(3) * d3f(x)
                                       * d3f(x))
                   / (-d0f(x) * d6f(x) - (long double)(2) * d1f(x) * d5f(x) + d2f(x) * d4f(x) + (long double)(2) *
                      d3f(x) * d3f(x));

        case 4: // После 5 применений Лопиталя
            return (long double)(2) * (-d1f(x) * d7f(x) - (long double)(6) * d2f(x) * d6f(x) - (long double)(15) *
                                       d3f(x) * d5f(x) - (long double)(10) * d4f(x) * d4f(x))
                   / (d0f(x) * d8f(x) + (long double)(4) * d1f(x) * d7f(x) + (long double)(4) * d2f(x) * d6f(x)
                      - (long double)(4) * d3f(x) * d5f(x) - (long double)(5) * d4f(x) * d4f(x));

        case 5:
            return (long double)(2) * (-d1f(x) * d9f(x) - (long double)(8) * d2f(x) * d8f(x)
                                       - (long double)(28) * d3f(x) * d7f(x) - (long double)(56) * d4f(x) * d6f(x)
                                       - (long double)(35) * d5f(x) * d5f(x))
                   / (d0f(x) * d10f(x) + (long double)(6) * d1f(x) * d9f(x) + (long double)(13) * d2f(x) * d8f(x)
                      + (long double)(8) * d3f(x) * d7f(x) - (long double)(14) * d4f(x) * d6f(x)
                      - (long double)(14) * d5f(x) * d5f(x));

        default:
            return (long double)(0); // или throw std::invalid_argument("m > 6 не поддерживается");
    }
}

Polynomial::Type Polynomial::multiplicity(const Complex& x, bool f) const {
    const auto d1f = derivative();
    const auto d2f = d1f.derivative();
    auto s         = d1f(x);
    s              *= s;
    return std::abs(Complex(1) / (Complex(1) - (*this)(x) * d2f(x) / s));
}

#include <iostream>

int Polynomial::computeMultiplicity(const Type& x) const {
    /*if (std::abs((*this)(x) / x) > 1e-3)
        return 1;*/

    auto num = derivative();
    num      *= num;
    auto den = num - *this * derivative(2);

    for (int m = 1; m < n; ++m) {
        std::cout << "m(" << m << ") = " << num(x) / den(x) << '\n';
        if (const auto res = num(x) / den(x); res > 0)
            if (std::abs(res - m) < 5e-2)
                return m;

        num.makeDerivative(2);
        den.makeDerivative(2);
    }

    return 1;
}

int Polynomial::computeMultiplicity(const Complex& x) const {
    if (std::abs((*this)(x) / x) > 1e-3)
        return 1;

    auto num = derivative();
    num      *= num;
    auto den = num - *this * derivative(2);

    for (int m = 1; m < n; ++m) {
        if (const auto res = num(x) / den(x); res.real() > 0)
            if (std::abs(std::abs(res) - m) < 5e-2)
                return m;

        num.makeDerivative(2);
        den.makeDerivative(2);
    }

    return 1;
}

int Polynomial::compute_multiplicity(const Type& x, double tol) const {
    const Type px = (*this)(x);
    //if (std::abs(px) > tol) return 0;

    const int n  = degree();
    double scale = 0.0;
    /*for (int k = 0; k <= n; ++k) {
        scale = std::max(scale, std::abs(derivative(k)(x)));
    }*/
    if (scale == 0.0)
        return n + 1; // нулевой полином

    int m = 0;
    for (int k = 0; k <= n; ++k) {
        if (std::abs(derivative(k)(x)) > tol * scale)
            break;
        m = k + 1;
    }
    return m;
}
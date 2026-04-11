//
// Created by Vadim on 21.08.2025.
//

#include "../../include/numina/polynomial.h"

Polynomial::Type Polynomial::multiplicity(const Type& x) const {
    const auto d1f = derivative();
    const auto d2f = d1f.derivative();
    auto s = d1f(x); s *= s;
    return std::abs(s / (s - (*this)(x) * d2f(x)));
}

Polynomial::Type Polynomial::multiplicity(const Complex& x) const {
    const auto d1f = derivative();
    const auto d2f = d1f.derivative();
    auto s = d1f(x); s *= s;
    return std::abs(s / (s - (*this)(x) * d2f(x)));
}

int Polynomial::computeMultiplicity(const Type& x) const {
    if (std::abs((*this)(x) / x) > 1e-3)
        return 1;

    auto num = derivative(); num *= num;
    auto den = num - *this * derivative(2);

    for (int m = 1; m < n; ++m) {
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

    auto num = derivative(); num *= num;
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
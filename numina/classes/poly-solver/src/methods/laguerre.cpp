#include "numina/poly-solver.h"

// Created by Vadim on 21.08.2025.
namespace numina {
PolySolver::Complex PolySolver::polish_explicit_laguerre(Complex x) const noexcept {
    Complex fx = f(x);
    for (std::size_t i = 0; i < LAGUERRE_MAX_ITER; ++i) {
        if (std::abs(fx) < E1 * std::abs(x))
            return x;

        const Complex d1fx = d1(x);
        if (std::abs(d1fx) < 1e-20 * std::abs(fx)) {
            auto h  = static_cast<Type>(degree) * fx;
            auto f1 = f(x + h);
            auto f2 = f(x - h);
            x       += std::abs(f1) < std::abs(f2) ? h : -h;
            fx      = f(x);
            continue;
        }
        const Complex b            = static_cast<Type>(degree - 1) * d1fx;
        const Complex discriminant = std::sqrt(b * b - static_cast<Type>(degree * (degree - 1)) * fx * d2(x));
        const Complex a1           = d1fx - discriminant;
        const Complex a2           = d1fx + discriminant;

        x  -= static_cast<Type>(degree) * fx / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        fx = f(x);
    }

    return x;
}

PolySolver::Complex PolySolver::polish_implicit_laguerre(Complex x) const noexcept {
    Complex fx = f(x);
    for (std::size_t i = 0; i < LAGUERRE_MAX_ITER; ++i) {
        if (std::abs(fx) < E1 * std::abs(x))
            return x;
        const Complex d1fx = d1(x);
        if (std::abs(d1fx) < 1e-20 * std::abs(fx)) {
            auto h  = static_cast<Type>(degree) * fx;
            auto f1 = f(x + h);
            auto f2 = f(x - h);
            x       += std::abs(f1) < std::abs(f2) ? h : -h;
            fx      = f(x);
            continue;
        }

        Complex S(0.0), T(0.0);
        for (const auto& [mu, r] : found) {
            const Complex den = x - r;
            S                 += static_cast<Type>(mu) / den;
            T                 += static_cast<Type>(mu) / (den * den);
        }

        const Complex G0 = d1fx / fx;
        const Complex H0 = G0 * G0 - d2(x) / fx;
        const Complex G  = G0 - S;
        const Complex H  = H0 - T;
        const Complex sd = std::sqrt(static_cast<Type>(m_eff * (m_eff - 1)) * H - G * G);
        const Complex a1 = G + sd;
        const Complex a2 = G - sd;

        x  -= static_cast<Type>(m_eff) / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        fx = f(x);
    }

    return x;
}
}
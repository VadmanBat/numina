#include "numina/poly-solver.h"

// Created by Vadim on 22.04.2026.
namespace numina {
PolySolver::Complex PolySolver::polish_explicit_newton(Complex x, const std::size_t m) const {
    const auto& g  = df[m - 1];
    const auto& dg = df[m];

    Complex gx = g(x);
    for (std::size_t i = 0; i < NEWTON_MAX_ITER; ++i) {
        Complex step = gx / dg(x);
        x            -= step;
        gx           = g(x);

        if (std::abs(step) < E1 * (Type(1.0) + std::abs(x)))
            break;
    }

    return x;
}

PolySolver::Complex PolySolver::polish_implicit_newton(Complex x, const std::size_t m) {
    const auto& g  = df[m - 1];
    const auto& dg = df[m];

    std::vector<std::pair<std::size_t, Complex>> roots(found.lower_bound(m), found.end());

    Complex gx = g(x);
    for (std::size_t i = 0; i < NEWTON_MAX_ITER; ++i) {
        Complex dg_val = dg(x);

        Complex S(0.0);
        for (const auto& [mu, r] : roots)
            S += static_cast<Type>(mu - m + 1) / (x - r);

        Complex denom = dg_val - gx * S;
        if (std::abs(denom) < E2)
            break;

        Complex step = gx / denom;
        x            -= step;
        gx           = g(x);

        if (std::abs(step) < E1 * (Type(1.0) + std::abs(x)))
            break;
    }

    return x;
}
}
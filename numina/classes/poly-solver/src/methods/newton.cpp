#include "numina/poly-solver.h"

// Created by Vadim on 22.04.2026.
namespace numina {
PolySolver::Complex PolySolver::polish_explicit_newton(Complex x, const std::size_t m) const {
    const auto& g  = df[m - 1];
    const auto& dg = df[m];

    Complex gx    = g(x);
    Type prev_abs = std::abs(gx);
    if (prev_abs == 0)
        return x;

    for (std::size_t i = 0; i < NEWTON_MAX_ITER; ++i) {
        Complex x_new      = x - gx / dg(x);
        Complex gx_new     = g(x_new);
        const Type new_abs = std::abs(gx_new);

        if (new_abs >= prev_abs)
            return x;

        x        = x_new;
        gx       = gx_new;
        prev_abs = new_abs;

        if (new_abs == 0)
            return x;
    }

    return x;
}

PolySolver::Complex PolySolver::polish_implicit_newton(Complex x, const std::size_t m) const {
    const auto& g  = df[m - 1];
    const auto& dg = df[m];

    std::vector<std::pair<std::size_t, Complex>> roots(found.lower_bound(m), found.end());

    Complex gx    = g(x);
    Type prev_abs = std::abs(gx);
    if (prev_abs == 0)
        return x;

    for (std::size_t i = 0; i < NEWTON_MAX_ITER; ++i) {
        Complex S(0.0);
        for (const auto& [mu, r] : roots)
            S += static_cast<Type>(mu - m + 1) / (x - r);

        Complex x_new      = x - gx / (dg(x) - gx * S);
        Complex gx_new     = g(x_new);
        const Type new_abs = std::abs(gx_new);

        if (new_abs >= prev_abs)
            return x;

        x        = x_new;
        gx       = gx_new;
        prev_abs = new_abs;

        if (new_abs == 0)
            return x;
    }

    return x;
}
}
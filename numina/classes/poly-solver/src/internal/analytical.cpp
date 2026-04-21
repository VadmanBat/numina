#include "numina/poly-solver.h"
// Created by Vadim on 19.04.2026.
namespace numina {
void PolySolver::solve_general_case() {
    std::size_t m_eff = degree;
    while (m_eff != 0) {
        auto x = laguerre(Complex(-1.0, -1.0));
        auto m = df[0].computeMultiplicity(x);
        if (m > 2)
            for (std::size_t i = df.size() - 1; i < m; ++i)
                df.emplace_back(df[i].derivative());
        x = newton(x, m);

        found.emplace_back(x, m);
        m_eff -= m;

        if (std::abs(x.imag()) > E1 * std::abs(x)) {
            Complex conj_x = std::conj(x);
            found.emplace_back(conj_x, m);
            m_eff -= m;

            result.second.emplace_back(x, m);
            result.second.emplace_back(conj_x, m);
        }
        else
            result.first.emplace_back(x.real(), m);
    }
}

void PolySolver::solve_quadratic_case() {
#define a coeffs[0]
#define b coeffs[1]
#define c coeffs[2]

    if (const auto D = b * b - 4 * a * c; D == 0)
        result.first.emplace_back(-0.5 * b / a, 2);
    else if (D > 0) {
        const auto q = b >= 0
                           ? -0.5 * (b + std::sqrt(D))
                           : -0.5 * (b - std::sqrt(D));
        result.first.emplace_back(q / a, 1);
        result.first.emplace_back(c / q, 1);
    }
    else {
        const auto q = b >= 0
                           ? -0.5 * (b + std::sqrt(Complex(D)))
                           : -0.5 * (b - std::sqrt(Complex(D)));
        result.second.emplace_back(q / a, 1);
        result.second.emplace_back(c / q, 1);
    }

#undef a
#undef b
#undef c
}

void PolySolver::solve_cases() {
    std::size_t zero_mult = 0;
    while (zero_mult < degree && coeffs[degree - zero_mult] == 0)
        ++zero_mult;

    if (zero_mult != 0) {
        result.first.emplace_back(0, zero_mult);
        coeffs.resize(coeffs.size() - zero_mult);
        degree -= zero_mult;
    }

    switch (degree) {
        case 0:
            break;
        case 1:
            result.first.emplace_back(-coeffs[1] / coeffs[0], 1);
            break;
        case 2:
            solve_quadratic_case();
            break;
        default:
            //solve_general_case();
            solve_with_implicit_deflation();
    }

    clear();
}
}
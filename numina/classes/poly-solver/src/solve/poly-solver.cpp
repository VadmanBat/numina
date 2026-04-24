#include "numina/poly-solver.h"
// Created by Vadim on 19.04.2026.
namespace numina {
void PolySolver::solve_explicit_general_case() noexcept {
    while (degree != 0) {
        auto x = polish_explicit_laguerre(0);
        auto m = df[0].computeMultiplicity(x);
        if (m > 2)
            for (std::size_t i = df.size() - 1; i < m; ++i)
                df.emplace_back(df[i].derivative());
        x = polish_explicit_newton(x, m);

        if (std::abs(x.imag()) > E1 * std::abs(x)) {
            Complex conj_x = std::conj(x);

            answer.second.emplace_back(x, m);
            answer.second.emplace_back(conj_x, m);

            deflate_conj(x, m);
        }
        else {
            answer.first.emplace_back(x.real(), m);
            deflate(x.real(), m);
        }
    }
}

void PolySolver::solve_implicit_general_case() noexcept {
    m_eff = degree;
    while (m_eff != 0) {
        auto x = polish_implicit_laguerre(0);
        auto m = df[0].computeMultiplicity(x);
        if (m > 2)
            for (std::size_t i = df.size() - 1; i < m; ++i)
                df.emplace_back(df[i].derivative());
        x = polish_explicit_newton(x, m);

        found.emplace(m, x);
        m_eff -= m;

        if (std::abs(x.imag()) > E1 * std::abs(x)) {
            Complex conj_x = std::conj(x);
            found.emplace(m, conj_x);
            m_eff -= m;

            answer.second.emplace_back(x, m);
            answer.second.emplace_back(conj_x, m);
        }
        else
            answer.first.emplace_back(x.real(), m);
    }
}

void PolySolver::solve_quadratic_case() noexcept {
#define a coeffs[0]
#define b coeffs[1]
#define c coeffs[2]

    if (const auto D = b * b - 4 * a * c; D == 0)
        answer.first.emplace_back(-0.5 * b / a, 2);
    else if (D > 0) {
        const auto q = b >= 0
                           ? -0.5 * (b + std::sqrt(D))
                           : -0.5 * (b - std::sqrt(D));
        answer.first.emplace_back(q / a, 1);
        answer.first.emplace_back(c / q, 1);
    }
    else {
        const auto q = b >= 0
                           ? -0.5 * (b + std::sqrt(Complex(D)))
                           : -0.5 * (b - std::sqrt(Complex(D)));
        const Complex r1 = q / a;
        answer.second.emplace_back(r1, 1);
        answer.second.emplace_back(std::conj(r1), 1);
    }

#undef a
#undef b
#undef c
}

void PolySolver::solve_cases() noexcept {
    std::size_t zero_mult = 0;
    while (zero_mult < degree && coeffs[degree - zero_mult] == 0)
        ++zero_mult;

    if (zero_mult != 0) {
        answer.first.emplace_back(0, zero_mult);
        coeffs.resize(coeffs.size() - zero_mult);
        degree -= zero_mult;
    }

    switch (degree) {
        case 0:
            break;
        case 1:
            answer.first.emplace_back(-coeffs[1] / coeffs[0], 1);
            break;
        case 2:
            solve_quadratic_case();
            break;
        default:
            solve_implicit_general_case();
    }

    clear();
}
}
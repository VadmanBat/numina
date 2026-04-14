#include "numina/quadratic-deriv-poly.h"
// Created by Vadim on 14.04.2026.
#include <algorithm>

namespace numina {
QuadraticDerivPoly::QuadraticDerivPoly(std::initializer_list<std::tuple<int, int, ll>> monomials) {
    for (const auto& [i, j, c] : monomials)
        addTerm(i, j, c);
}

void QuadraticDerivPoly::addTerm(int i, int j, const ll coeff) {
    if (i > j)
        std::swap(i, j);
    if (coeff != 0)
        terms[{i, j}] += coeff;
    else
        terms.erase({i, j});
}

const QuadraticDerivPoly::Terms& QuadraticDerivPoly::getTerms() const noexcept {
    return terms;
}

void QuadraticDerivPoly::differentiate() {
    Terms nt;
    for (const auto& [pr, c] : terms) {
        const auto [a, b]         = pr;
        nt[std::minmax(a + 1, b)] += c; // p^{(a+1)} p^{(b)}
        nt[{a, b + 1}]            += c; // p^{(a)} p^{(b+1)}
    }
    terms = std::move(nt);
} // d/dt (p^{(a)} p^{(b)}) = p^{(a+1)} p^{(b)} + p^{(a)} p^{(b+1)}

void QuadraticDerivPoly::doubleDifferentiate() {
    Terms nt;
    for (const auto& [pr, c] : terms) {
        const auto [a, b]         = pr;
        nt[std::minmax(a + 2, b)] += c;     // p^(a+2)p^b
        nt[{a + 1, b + 1}]        += 2 * c; // 2 * p^(a+1)p^(b+1)
        nt[{a, b + 2}]            += c;     // p^a p^(b+2)
    }
    terms = std::move(nt);
} // d²/dt² (p^{(a)} p^{(b)}) = p^{(a+2)} p^{(b)} + 2 * p^(a+1)p^(b+1) + p^{(a)} p^{(b+2)}
}
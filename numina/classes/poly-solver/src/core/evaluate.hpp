#pragma once
// Created by Vadim on 13.07.2025.
namespace numina {
inline PolySolver::Complex PolySolver::f(const Complex& x) const noexcept {
    Complex result = c[0];
    for (std::size_t i = 1; i <= degree; ++i)
        (result *= x) += c[i];
    return result;
}

inline PolySolver::Complex PolySolver::d1(const Complex& x) const noexcept {
    Complex result = c_d1[0];
    for (std::size_t i = 1; i < degree; ++i)
        (result *= x) += c_d1[i];
    return result;
}

inline PolySolver::Complex PolySolver::d2(const Complex& x) const noexcept {
    const auto n   = degree - 1;
    Complex result = c_d2[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c_d2[i];
    return result;
}
}
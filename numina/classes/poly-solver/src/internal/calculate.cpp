#include "numina/poly-solver.h"
// Created by Vadim on 13.07.2025.
namespace numina {
PolySolver::Complex PolySolver::f(const Complex& x) const {
    Complex result = c[0];
    for (std::size_t i = 1; i <= degree; ++i)
        (result *= x) += c[i];
    return result;
}

PolySolver::Complex PolySolver::d1(const Complex& x) const {
    Complex result = c_d1[0];
    for (std::size_t i = 1; i < degree; ++i)
        (result *= x) += c_d1[i];
    return result;
}

PolySolver::Complex PolySolver::d2(const Complex& x) const {
    const auto n = degree - 1;
    Complex result = c_d2[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c_d2[i];
    return result;
}

int PolySolver::compute_multiplicity(const Complex& x) {
    int m = 1;
    const auto n = static_cast<int>(nums.size());
    for (; m < n; ++m) {
        const auto& num = nums[m - 1];
        const auto& den = dens[m - 1];

        if (const auto res = num(x) / den(x); res.real() > 0)
            if (std::abs(std::abs(res) - m) < 5e-2)
                return m;
    }

    for (; m < degree; ++m) {
        const auto& num = nums[m - 1];
        const auto& den = dens[m - 1];

        if (const auto res = num(x) / den(x); res.real() > 0)
            if (std::abs(std::abs(res) - m) < 5e-2)
                return m;

        nums.emplace_back(num.derivative(2));
        dens.emplace_back(den.derivative(2));
    }

    return 1;
}
}
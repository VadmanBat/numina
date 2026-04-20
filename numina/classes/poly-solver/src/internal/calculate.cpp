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

std::size_t PolySolver::compute_multiplicity(const Complex& x) {
    std::vector<Complex> dfx_1;
    std::vector<LongComplex> dfx_2;
    dfx_1.reserve(degree);
    dfx_2.reserve(degree);

    auto deriv = *this;
    for (int i = 0; i < 3; ++i) {
        dfx_1.push_back(df[i](static_cast<Complex>(x)));
        dfx_2.push_back(df[i](static_cast<LongComplex>(x)));
    }

    int prev_m = -1;
    for (std::size_t m = 1; m <= degree; ++m) {
        const Type m1 = static_cast<Type>(m) - 1 + (dfx_1[m] * dfx_1[m] / (
                               dfx_1[m] * dfx_1[m] - dfx_1[m - 1] * dfx_1[m + 1])).real();
        const LongType m2 = static_cast<LongType>(m) - 1 + (dfx_2[m] * dfx_2[m] / (
                                   dfx_2[m] * dfx_2[m] - dfx_2[m - 1] * dfx_2[m + 1])).real();

        if (std::abs(m1 - m2) < 1e-6)
            if (m == std::round(m2))
                return m;
        if (const auto rounded = static_cast<int>(std::round(m1)); rounded == std::round(m2)) {
            if (prev_m >= static_cast<int>(m) && prev_m == rounded)
                return prev_m;
            prev_m = rounded;
        }

        if (df.size() == m + 2)
            df.emplace_back(df[m + 1].derivative());
        dfx_1.push_back(df[m + 2](static_cast<Complex>(x)));
        dfx_2.push_back(df[m + 2](static_cast<LongComplex>(x)));
    }

    return 1;
}
}
#include "numina/poly-solver.h"
// Created by Vadim on 18.08.2025.
namespace numina {
void PolySolver::trim_leading_zeros() noexcept {
    if (coeffs.empty()) {
        coeffs.assign(1, 0);
        return;
    }

    if (coeffs[0] != 0)
        return;

    std::size_t first_non_zero = 1;
    const auto n               = coeffs.size();
    while (first_non_zero < n && coeffs[first_non_zero] == 0)
        ++first_non_zero;

    if (first_non_zero == n) {
        coeffs.assign(1, 0);
        return;
    }

    std::move(coeffs.begin() + static_cast<long long>(first_non_zero),
              coeffs.end(),
              coeffs.begin());
    coeffs.resize(n - first_non_zero);
}

void PolySolver::prepare() {
    degree = coeffs.size() - 1;

    Type lead = std::abs(coeffs[0]);
    for (auto& c : coeffs)
        c /= lead;

    coeffs_d1.resize(degree);
    coeffs_d2.resize(degree - 1);

    c    = coeffs.data();
    c_d1 = coeffs_d1.data();
    c_d2 = coeffs_d2.data();

    compute_derivative();

    df.reserve(degree + 1);
    df.emplace_back(coeffs);
    df.emplace_back(coeffs_d1);
    df.emplace_back(coeffs_d2);

    /*nums.reserve(degree);
    dens.reserve(degree);
    auto num = df[1] * df[1];
    auto den = num - df[0] * df[2];
    nums.emplace_back(std::move(num));
    dens.emplace_back(std::move(den));*/
}

void PolySolver::prepare(const std::vector<Type>& coefficients) {
    coeffs = coefficients;
    trim_leading_zeros();
    prepare();
}

void PolySolver::prepare(std::vector<Type>&& coefficients) {
    coeffs = std::move(coefficients);
    trim_leading_zeros();
    prepare();
}

void PolySolver::prepare(const Polynomial& poly) {
    coeffs = poly.vector();
    prepare();
}

void PolySolver::prepare(Polynomial&& poly) {
    coeffs = std::move(poly).extractVector();
    prepare();
}

void PolySolver::clear() {
    df.clear();
    /*nums.clear();
    dens.clear();*/
}
}
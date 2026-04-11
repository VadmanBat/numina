#include "numina/poly-solver.h"
// Created by Vadim on 18.08.2025.
namespace numina {
void PolySolver::setup_state(const std::vector <Type>& coefficients) {
    coeffs = coefficients;
    degree = coeffs.size() - 1;

    coeffs_d1.resize(degree);
    coeffs_d2.resize(degree - 1);

    c = coeffs.data();
    c_d1 = coeffs_d1.data();
    c_d2 = coeffs_d2.data();

    compute_derivative();

    dfx.reserve(degree + 1);
    dfx.emplace_back(coeffs);
    dfx.emplace_back(coeffs_d1);
    dfx.emplace_back(coeffs_d2);

    nums.reserve(degree);
    dens.reserve(degree);
    auto num = dfx[1] * dfx[1];
    auto den = num - dfx[0] * dfx[2];
    nums.emplace_back(std::move(num));
    dens.emplace_back(std::move(den));
}

void PolySolver::setup_state(std::vector <Type>&& coefficients) {
    coeffs = std::move(coefficients);
    degree = coeffs.size() - 1;

    coeffs_d1.resize(degree);
    coeffs_d2.resize(degree - 1);

    c = coeffs.data();
    c_d1 = coeffs_d1.data();
    c_d2 = coeffs_d2.data();

    compute_derivative();

    dfx.reserve(degree + 1);
    dfx.emplace_back(coeffs);
    dfx.emplace_back(coeffs_d1);
    dfx.emplace_back(coeffs_d2);

    nums.reserve(degree);
    dens.reserve(degree);
    auto num = dfx[1] * dfx[1];
    auto den = num - dfx[0] * dfx[2];
    nums.emplace_back(std::move(num));
    dens.emplace_back(std::move(den));
}

void PolySolver::setup_state(const Polynomial& poly) {
    coeffs = poly.vector();
    degree = coeffs.size() - 1;

    coeffs_d1.resize(degree);
    coeffs_d2.resize(degree - 1);

    c = coeffs.data();
    c_d1 = coeffs_d1.data();
    c_d2 = coeffs_d2.data();

    compute_derivative();

    dfx.reserve(degree + 1);
    dfx.emplace_back(coeffs);
    dfx.emplace_back(coeffs_d1);
    dfx.emplace_back(coeffs_d2);

    nums.reserve(degree);
    dens.reserve(degree);
    auto num = dfx[1] * dfx[1];
    auto den = num - dfx[0] * dfx[2];
    nums.emplace_back(std::move(num));
    dens.emplace_back(std::move(den));
}

void PolySolver::setup_state(Polynomial&& poly) {
    coeffs = poly.vector();
    degree = coeffs.size() - 1;

    coeffs_d1.resize(degree);
    coeffs_d2.resize(degree - 1);

    c = coeffs.data();
    c_d1 = coeffs_d1.data();
    c_d2 = coeffs_d2.data();

    compute_derivative();

    dfx.reserve(degree + 1);
    dfx.emplace_back(std::move(poly));
    dfx.emplace_back(coeffs_d1);
    dfx.emplace_back(coeffs_d2);

    nums.reserve(degree);
    dens.reserve(degree);
    auto num = dfx[1] * dfx[1];
    auto den = num - dfx[0] * dfx[2];
    nums.emplace_back(std::move(num));
    dens.emplace_back(std::move(den));
}

void PolySolver::clear_state() {
    dfx.clear();
    nums.clear();
    dens.clear();
}
    }
#pragma once
// Created by Vadim on 12.07.2025.
#include "numina/polynomial.h"

#include <vector>
#include <map>
#include <complex>
#include <limits>
#include <utility>

namespace numina {
class PolySolver {
public:
    using Type        = double;
    using Complex     = std::complex<Type>;
    using LongType    = long double;
    using LongComplex = std::complex<LongType>;
    using Roots       = std::pair<
        std::vector<std::pair<Type, int>>,
        std::vector<std::pair<Complex, int>>
    >;

    PolySolver() = default;

    std::vector<Complex> solve(const std::vector<Type>& coefficients);
    std::vector<Complex> solve(std::vector<Type>&& coefficients);
    std::vector<Complex> solve(const Polynomial& poly);
    std::vector<Complex> solve(Polynomial&& poly);

    Roots solveWithMultiplicities(const std::vector<Type>& coefficients);
    Roots solveWithMultiplicities(std::vector<Type>&& coefficients);
    Roots solveWithMultiplicities(const Polynomial& poly);
    Roots solveWithMultiplicities(Polynomial&& poly);

private:
    inline static const Type E1                    = std::sqrt(std::numeric_limits<Type>::epsilon());
    inline static const Type E2                    = std::sqrt(E1);
    static constexpr std::size_t LAGUERRE_MAX_ITER = 30;
    static constexpr std::size_t NEWTON_MAX_ITER   = 30;

    std::size_t degree = 0;

    std::vector<Type> coeffs;
    std::vector<Type> coeffs_d1;
    std::vector<Type> coeffs_d2;

    Type* c    = nullptr;
    Type* c_d1 = nullptr;
    Type* c_d2 = nullptr;

    std::vector<Polynomial> df;
    std::multimap<std::size_t, Complex> found;

    Roots answer;

    void trim_leading_zeros() noexcept;
    void prepare();
    void prepare(const std::vector<Type>& coefficients);
    void prepare(std::vector<Type>&& coefficients);
    void prepare(const Polynomial& poly);
    void prepare(Polynomial&& poly);
    void clear();

    void compute_derivative() const;

    [[nodiscard]] [[gnu::always_inline]] Complex f(const Complex& x) const;
    [[nodiscard]] [[gnu::always_inline]] Complex d1(const Complex& x) const;
    [[nodiscard]] [[gnu::always_inline]] Complex d2(const Complex& x) const;

    void deflate(const Type& root, std::size_t m = 1);
    void deflate_conj(const Complex& root, std::size_t m = 1);

    [[nodiscard]] Complex polish_explicit_laguerre(Complex x) const;
    [[nodiscard]] Complex polish_implicit_laguerre(Complex x) const;
    [[nodiscard]] Complex polish_explicit_newton(Complex x, std::size_t m) const;
    [[nodiscard]] Complex polish_implicit_newton(Complex x, std::size_t m) const;

    void solve_explicit_general_case();
    void solve_implicit_general_case();
    void solve_quadratic_case();
    void solve_cases();

    std::vector<Complex> get_vector();
};

inline PolySolver::Complex PolySolver::f(const Complex& x) const {
    Complex result = c[0];
    for (std::size_t i = 1; i <= degree; ++i)
        (result *= x) += c[i];
    return result;
}

inline PolySolver::Complex PolySolver::d1(const Complex& x) const {
    Complex result = c_d1[0];
    for (std::size_t i = 1; i < degree; ++i)
        (result *= x) += c_d1[i];
    return result;
}

inline PolySolver::Complex PolySolver::d2(const PolySolver::Complex& x) const {
    const auto n   = degree - 1;
    Complex result = c_d2[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= x) += c_d2[i];
    return result;
}
}

/*
 *
 * расчёты
 *
Перспективы повышения точности:
    - вычисление значений полинома в точках: f(x);
    - вычисление производных: f'(x), f''(x), ...;
    - дефляция основного полинома: f(x) /= (x - r)^m
*/
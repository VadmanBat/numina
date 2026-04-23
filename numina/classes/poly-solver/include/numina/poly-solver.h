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

    Roots solveWithImplicitDeflation(const std::vector<Type>& coefficients); // test-legacy
    Roots solveWithImplicitDeflation(std::vector<Type>&& coefficients); // test-legacy
    Roots solveWithImplicitDeflation(const Polynomial& poly); // test-legacy
    Roots solveWithImplicitDeflation(Polynomial&& poly); // test-legacy

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

    Roots result;

    void compute_derivative() const;

    void trim_leading_zeros() noexcept;
    void prepare();
    void prepare(const std::vector<Type>& coefficients);
    void prepare(std::vector<Type>&& coefficients);
    void prepare(const Polynomial& poly);
    void prepare(Polynomial&& poly);
    void clear();

    [[nodiscard]] Complex f(const Complex& x) const;
    [[nodiscard]] Complex d1(const Complex& x) const;
    [[nodiscard]] Complex d2(const Complex& x) const;

    std::size_t compute_multiplicity(const Complex& x); // test-legacy
    void deflate(const Type& root, int m = 1);
    void deflate_conj(const Complex& root, int m = 1);

    std::vector<Complex> solve(); // test-legacy
    void solve_with_multiplicities(); // test-legacy
    void solve_with_implicit_deflation(); // test-legacy

    Complex polish_implicit_laguerre(Complex x) {}
    Complex polish_explicit_laguerre(Complex x) {}
    Complex polish_implicit_newton(Complex x, std::size_t m) const;
    Complex polish_explicit_newton(Complex x, std::size_t m) const;

    void solve_general_case();
    void solve_quadratic_case();
    void solve_cases();
};
}

/*
 * расчёт кратности!
 * логика расчёта!
 * расчёты
 *
Перспективы повышения точности:
    - вычисление значений полинома в точках: f(x);
    - вычисление производных: f'(x), f''(x), ...;
    - стабильное нахождение кратности корня: m;
    - дефляция основного полинома: f(x) /= (x - r)^m
*/
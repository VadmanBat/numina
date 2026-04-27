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
    using Type       = double;
    using Complex    = std::complex<Type>;
    using MultiRoots = std::pair<
        std::vector<std::pair<Type, std::size_t>>,
        std::vector<std::pair<Complex, std::size_t>>
    >;

    enum class Method : bool {
        Explicit = false,
        Implicit = true
    };

    PolySolver() = default;

    std::vector<Complex> solve(const std::vector<Type>& coefficients, Method method = Method::Explicit);
    std::vector<Complex> solve(std::vector<Type>&& coefficients, Method method = Method::Explicit);
    std::vector<Complex> solve(const Polynomial& poly, Method method = Method::Explicit);
    std::vector<Complex> solve(Polynomial&& poly, Method method = Method::Explicit);

    MultiRoots multisolve(const std::vector<Type>& coefficients, Method method = Method::Explicit);
    MultiRoots multisolve(std::vector<Type>&& coefficients, Method method = Method::Explicit);
    MultiRoots multisolve(const Polynomial& poly, Method method = Method::Explicit);
    MultiRoots multisolve(Polynomial&& poly, Method method = Method::Explicit);

private:
    inline static const Type E1                    = std::sqrt(std::numeric_limits<Type>::epsilon());
    inline static const Type E2                    = std::sqrt(E1);
    static constexpr std::size_t LAGUERRE_MAX_ITER = 30;
    static constexpr std::size_t NEWTON_MAX_ITER   = 30;

    std::size_t degree = 0, m_eff = 0;

    std::vector<Type> coeffs;
    std::vector<Type> coeffs_d1;
    std::vector<Type> coeffs_d2;

    Type* c    = nullptr;
    Type* c_d1 = nullptr;
    Type* c_d2 = nullptr;

    std::vector<Polynomial> df;
    std::multimap<std::size_t, Complex> found;

    MultiRoots answer;

    void trim_leading_zeros() noexcept;
    void prepare();
    void prepare(const std::vector<Type>& coefficients);
    void prepare(std::vector<Type>&& coefficients);
    void prepare(const Polynomial& poly);
    void prepare(Polynomial&& poly);
    void clear() noexcept;

    void compute_derivative() const noexcept;

    [[nodiscard]] [[gnu::always_inline]] Complex f(const Complex& x) const noexcept;
    [[nodiscard]] [[gnu::always_inline]] Complex d1(const Complex& x) const noexcept;
    [[nodiscard]] [[gnu::always_inline]] Complex d2(const Complex& x) const noexcept;

    void deflate(const Type& root, std::size_t m = 1) noexcept;
    void deflate_conj(const Complex& root, std::size_t m = 1) noexcept;

    [[nodiscard]] Complex polish_explicit_laguerre(Complex x) const noexcept;
    [[nodiscard]] Complex polish_implicit_laguerre(Complex x) const noexcept;
    [[nodiscard]] Complex polish_explicit_newton(Complex x, std::size_t m) const noexcept;
    [[nodiscard]] Complex polish_implicit_newton(Complex x, std::size_t m) const noexcept;

    void solve_explicit_general_case() noexcept;
    void solve_implicit_general_case() noexcept;
    void solve_quadratic_case() noexcept;
    void solve_cases(Method method) noexcept;

    std::vector<Complex> extract_answer() noexcept;
};
}

#include "../../src/core/evaluate.hpp"

/*
Перспективы повышения точности:
    - теорема о рациональных корнях;
    - вычисление производных: f'(x), f''(x), ...;
    - вычисление значений полинома в точках (возможно, через CompEA): f(x), f'(x), f''(x), ...;
    - дефляция основного полинома: f(x) /= (x - r)^m;
    - определение комплексного корня через уменьшение мнимой части (на порядки) при уточнении Ньютоном;
*/
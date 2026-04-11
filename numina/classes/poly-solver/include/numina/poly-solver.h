//
// Created by Vadim on 12.07.2025.
//

#ifndef NUMINA_POLY_SOLVER_H
#define NUMINA_POLY_SOLVER_H

#include "numina/classes/polynomial/include/numina/polynomial.h"

class PolySolver {
private:
    using Type = double;
    using Complex = std::complex <Type>;

    inline static const Type e1 = std::sqrt(std::numeric_limits<Type>::epsilon());
    inline static const Type e2 = std::sqrt(e1);

    std::size_t degree = 0;

    std::vector <Type> coeffs;
    std::vector <Type> coeffs_d1;
    std::vector <Type> coeffs_d2;

    Type* c = nullptr;
    Type* c_d1 = nullptr;
    Type* c_d2 = nullptr;

    std::vector <Polynomial> dfx;
    std::vector <Polynomial> nums, dens;

    void compute_derivative();

    void setup_state(const std::vector <Type>& coefficients);
    void setup_state(std::vector <Type>&& coefficients);
    void setup_state(const Polynomial& poly);
    void setup_state(Polynomial&& poly);
    void clear_state();

    Complex f(const Complex& x);
    Complex d1(const Complex& x);
    Complex d2(const Complex& x);

    int compute_multiplicity(const Complex& x);

    void deflate(const Type& root, int m = 1);
    void deflate_conj(const Complex& root, int m = 1);

    std::vector <Complex> solve();

    using Roots = std::pair <std::vector <std::pair <Type, int>>, std::vector <std::pair <Complex, int>>>;
    Roots solve_with_multiplicities();

public:
    PolySolver() = default;

    std::vector <Complex> solve(const std::vector <Type>& coefficients);
    std::vector <Complex> solve(std::vector <Type>&& coefficients);
    std::vector <Complex> solve(const Polynomial& poly);
    std::vector <Complex> solve(Polynomial&& poly);

    Roots solveWithMultiplicities(const std::vector <Type>& coefficients);
    Roots solveWithMultiplicities(std::vector <Type>&& coefficients);
    Roots solveWithMultiplicities(const Polynomial& poly);
    Roots solveWithMultiplicities(Polynomial&& poly);
};

#endif //NUMINA_POLY_SOLVER_H

/*
Перспективы повышения точности:
    - сначала искать реальные корни, потом мнимые;
    - вычисление значений полинома в точках: f(x);
    - вычисление производных: f'(x), f''(x), ...;
    - стабильное нахождение кратности корня: m;
    - дефляция основного полинома: f(x) /= (x - r)^m
*/
#include "numina/poly-solver.h"
// Created by Vadim on 13.07.2025.
#include <iostream>
namespace numina {
std::vector <PolySolver::Complex> PolySolver::solve(const std::vector <Type>& coefficients) {
    prepare(coefficients);
    solve_cases();
    return get_vector();
}

std::vector <PolySolver::Complex> PolySolver::solve(std::vector <Type>&& coefficients) {
    prepare(std::move(coefficients));
    solve_cases();
    return get_vector();
}

std::vector <PolySolver::Complex> PolySolver::solve(const Polynomial& poly) {
    prepare(poly);
    solve_cases();
    return get_vector();
}

std::vector <PolySolver::Complex> PolySolver::solve(Polynomial&& poly) {
    prepare(std::move(poly));
    solve_cases();
    return get_vector();
}

PolySolver::Roots PolySolver::solveWithMultiplicities(const std::vector <Type>& coefficients) {
    prepare(coefficients);
    solve_cases();
    return std::move(answer);
}

PolySolver::Roots PolySolver::solveWithMultiplicities(std::vector <Type>&& coefficients) {
    prepare(std::move(coefficients));
    solve_cases();
    return std::move(answer);
}

PolySolver::Roots PolySolver::solveWithMultiplicities(const Polynomial& poly) {
    prepare(poly);
    solve_cases();
    return std::move(answer);
}

PolySolver::Roots PolySolver::solveWithMultiplicities(Polynomial&& poly) {
    prepare(std::move(poly));
    solve_cases();
    return std::move(answer);
}
}
//
// Created by Vadim on 13.07.2025.
//

#include "../include/numina/poly-solver.h"

std::vector <PolySolver::Complex> PolySolver::solve(const std::vector <Type>& coefficients) {
    setup_state(coefficients);
    return solve();
}

std::vector <PolySolver::Complex> PolySolver::solve(std::vector <Type>&& coefficients) {
    setup_state(std::move(coefficients));
    return solve();
}

std::vector <PolySolver::Complex> PolySolver::solve(const Polynomial& poly) {
    setup_state(poly);
    return solve();
}

std::vector <PolySolver::Complex> PolySolver::solve(Polynomial&& poly) {
    setup_state(std::move(poly));
    return solve();
}

PolySolver::Roots PolySolver::solveWithMultiplicities(const std::vector <Type>& coefficients) {
    setup_state(coefficients);
    return solve_with_multiplicities();
}

PolySolver::Roots PolySolver::solveWithMultiplicities(std::vector <Type>&& coefficients) {
    setup_state(std::move(coefficients));
    return solve_with_multiplicities();
}

PolySolver::Roots PolySolver::solveWithMultiplicities(const Polynomial& poly) {
    setup_state(poly);
    return solve_with_multiplicities();
}

PolySolver::Roots PolySolver::solveWithMultiplicities(Polynomial&& poly) {
    setup_state(std::move(poly));
    return solve_with_multiplicities();
}

#include "numina/poly-solver.h"
// Created by Vadim on 13.07.2025.
namespace numina {
std::vector <PolySolver::Complex> PolySolver::solve(const std::vector <Type>& coefficients, const Method method) {
    prepare(coefficients);
    solve_cases(method);
    return extract_answer();
}

std::vector <PolySolver::Complex> PolySolver::solve(std::vector <Type>&& coefficients, const Method method) {
    prepare(std::move(coefficients));
    solve_cases(method);
    return extract_answer();
}

std::vector <PolySolver::Complex> PolySolver::solve(const Polynomial& poly, const Method method) {
    prepare(poly);
    solve_cases(method);
    return extract_answer();
}

std::vector <PolySolver::Complex> PolySolver::solve(Polynomial&& poly, const Method method) {
    prepare(std::move(poly));
    solve_cases(method);
    return extract_answer();
}

PolySolver::MultiRoots PolySolver::multisolve(const std::vector <Type>& coefficients, const Method method) {
    prepare(coefficients);
    solve_cases(method);
    return std::move(answer);
}

PolySolver::MultiRoots PolySolver::multisolve(std::vector <Type>&& coefficients, const Method method) {
    prepare(std::move(coefficients));
    solve_cases(method);
    return std::move(answer);
}

PolySolver::MultiRoots PolySolver::multisolve(const Polynomial& poly, const Method method) {
    prepare(poly);
    solve_cases(method);
    return std::move(answer);
}

PolySolver::MultiRoots PolySolver::multisolve(Polynomial&& poly, const Method method) {
    prepare(std::move(poly));
    solve_cases(method);
    return std::move(answer);
}
}
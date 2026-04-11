//
// Created by Vadim on 13.07.2025.
//

#include "../../include/numina/polynomial.h"

#include "../../../poly-solver/include/numina/poly-solver.h"

std::vector <Polynomial::Complex> Polynomial::computeRoots() const {
    numina::PolySolver solver;
    return solver.solve(*this);
}

Polynomial::Roots Polynomial::computeRootsWithMultiplicity() const {
    numina::PolySolver solver;
    return solver.solveWithMultiplicities(*this);
}

bool Polynomial::isRoot(const Type& root, const Type& tolerance) const {
    return std::abs((*this)(root)) < tolerance;
}

bool Polynomial::isRoot(const Complex& root, const Type& tolerance) const {
    return std::abs((*this)(root)) < tolerance;
}

#define MAX_ITER 50

Polynomial::Type Polynomial::refineRoot(Type root, const Type& tolerance) const {
    const auto dx = derivative();
    for (int i = 0; i < MAX_ITER; ++i) {
        if (isRoot(root, tolerance))
            return root;
        root -= (*this)(root) / dx(root);
    }
    return root;
}

Polynomial::Complex Polynomial::refineRoot(Complex root, const Type& tolerance) const {
    const auto dx = derivative();
    for (int i = 0; i < MAX_ITER; ++i) {
        if (isRoot(root, tolerance))
            return root;
        root -= (*this)(root) / dx(root);
    }
    return root;
}

#undef MAX_ITER
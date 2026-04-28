#include "numina/polynomial.h"
#include "numina/poly-solver.h"
// Created by Vadim on 13.07.2025.
std::vector<Polynomial::Complex> Polynomial::computeRoots() const {
    numina::PolySolver solver;
    return solver.solve(*this);
}

Polynomial::MultiRoots Polynomial::computeMultiRoots() const {
    numina::PolySolver solver;
    return solver.multisolve(*this);
}
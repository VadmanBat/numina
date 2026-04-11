//
// Created by Vadim on 11.07.2025.
//

#include "../include/numina/polynomial.h"

Polynomial Polynomial::compose(const Polynomial& other) const {
    if (isZero() || other.isZero())
        return {};

    Polynomial result = c[0];
    for (std::size_t i = 1; i < n; ++i)
        (result *= other) += c[i];
    return result;
}
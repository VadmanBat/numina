#include "numina/polynomial.h"
// Created by Vadim on 20.06.2025.
void Polynomial::remove_leading_zeros() {
    static const Type epsilon = std::sqrt(std::numeric_limits<Type>::epsilon());

    if (std::abs(c[0]) > epsilon)
        return;

    for (std::size_t i = 1; i < n; ++i)
        if (std::abs(c[i]) > epsilon) {
            coeffs.erase(coeffs.begin(), coeffs.begin() + static_cast<long long>(i));
            n -= i;
            return;
        }

    zeroing();
}

void Polynomial::zeroing() {
    coeffs.resize(1);
    c[0] = 0;
    n    = 1;
}
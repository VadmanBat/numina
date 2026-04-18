#include "numina/polynomial.h"
// Created by Vadim on 20.06.2025.
bool Polynomial::operator==(const Polynomial& other) const {
    if (n != other.n)
        return false;
    for (std::size_t i = 0; i < n; ++i)
        if (c[i] != other.c[i])
            return false;
    return true;
}

bool Polynomial::operator!=(const Polynomial& other) const {
    if (n != other.n)
        return true;
    for (std::size_t i = 0; i < n; ++i)
        if (c[i] != other.c[i])
            return true;
    return false;
}
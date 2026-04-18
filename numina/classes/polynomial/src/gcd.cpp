#include "numina/polynomial.h"
// Created by Vadim on 11.07.2025.
Polynomial gcd(Polynomial a, Polynomial b) {
    while (!b.isZero()) {
        Polynomial t = a % b;
        a            = std::move(b);
        b            = std::move(t);
    }
    return a;
}
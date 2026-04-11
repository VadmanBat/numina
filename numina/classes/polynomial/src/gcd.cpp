//
// Created by Vadim on 11.07.2025.
//

#include "../include/numina/polynomial.h"

Polynomial gcd(Polynomial a, Polynomial b) {
    while (!b.isZero()) {
        Polynomial t = a % b;
        a = std::move(b);
        b = std::move(t);
    }
    return a;
}
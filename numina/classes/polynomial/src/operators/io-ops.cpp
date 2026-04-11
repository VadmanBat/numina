//
// Created by Vadim on 20.06.2025.
//

#include "../../include/numina/polynomial.h"

std::istream& operator>>(std::istream& is, Polynomial& poly) {
    std::size_t n;
    is >> n;
    std::vector <Polynomial::Type> coefficients(n);
    auto c = coefficients.data();
    for (std::size_t i = 0; i < n; ++i)
        is >> c[i];
    poly = std::move(coefficients);
    return is;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& poly) {
    using Type = Polynomial::Type;

    if (poly.n == 1)
        return os << poly.c[poly.n - 1];

    if (poly.n == 2) {
        if (poly.c[0] != Type(0)) {
            if (poly.c[0] < 0)
                os << '-';
            if (const auto abs_value = std::abs(poly.c[0]); abs_value != Type(1))
                os << abs_value << " * ";
            os << 'x';
        }

        return poly.c[1] == Type(0) ? os :
               os << (poly.c[1] > Type(0) ? " + " : " - ")
                  << std::abs(poly.c[1]);
    }

    if (poly.c[0] < 0)
        os << '-';
    if (const auto abs_value = std::abs(poly.c[0]); abs_value != Type(1))
        os << abs_value << " * ";
    os << "x^" << poly.n - 1;

    for (std::size_t i = 1; i < poly.n - 2; ++i) {
        const auto coeff = poly.c[i];
        if (coeff == Type(0))
            continue;
        os << (coeff > Type(0) ? " + " : " - ");
        if (const auto abs_value = std::abs(coeff); abs_value != Type(1))
            os << abs_value << " * ";
        os << "x^" << poly.n - 1 - i;
    }

    if (poly.c[poly.n - 2] != Type(0)) {
        os << (poly.c[poly.n - 2] > Type(0) ? " + " : " - ");
        if (const auto abs_value = std::abs(poly.c[poly.n - 2]); abs_value != Type(1))
            os << abs_value << " * ";
        os << 'x';
    }

    return poly.c[poly.n - 1] == Type(0) ? os :
           os << (poly.c[poly.n - 1] > Type(0) ? " + " : " - ")
              << std::abs(poly.c[poly.n - 1]);
}
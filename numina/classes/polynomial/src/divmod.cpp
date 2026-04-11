//
// Created by Vadim on 19.06.2025.
//

#include "../include/numina/polynomial.h"

std::pair <Polynomial, Polynomial> Polynomial::divmod(const Polynomial& divisor) const {
    if (divisor.isZero())
        throw std::invalid_argument("Division by zero polynomial");

    if (n < divisor.n) // Когда степень делителя больше делимого
        return {Polynomial(), *this};

    const auto leading = divisor.c[0];
    const auto dn = n - divisor.n;

    auto remainder_coeffs = coeffs;
    std::vector <Type> quotient_coeffs(dn + 1, 0.0);
    const auto rem = remainder_coeffs.data();
    const auto quo = quotient_coeffs.data();

    for (std::size_t i = 0; i <= dn; ++i) {
        const auto factor = rem[i] / leading;
        quo[i] = factor;

        // Вычитаем текущий член частного, умноженный на делитель
        for (std::size_t j = 0; j < divisor.n; ++j)
            rem[i + j] -= factor * divisor.c[j];
    }

    // Формируем остаток (последние m коэффициентов)
    std::vector <Type> rem_coeffs;
    if (divisor.n > 1)
        rem_coeffs = std::vector<Type>(
                remainder_coeffs.begin() + int(dn + 1),
                remainder_coeffs.end()
        );

    return {std::move(quotient_coeffs), std::move(rem_coeffs)};
}

Polynomial Polynomial::operator/(const Polynomial& divisor) const {
    if (divisor.isZero())
        throw std::invalid_argument("Division by zero polynomial");

    if (n < divisor.n) // Когда степень делителя больше делимого
        return {};

    const auto leading = divisor.c[0];
    const auto dn = n - divisor.n;

    auto remainder_coeffs = coeffs;
    std::vector <Type> quotient_coeffs(dn + 1, 0.0);
    const auto rem = remainder_coeffs.data();
    const auto quo = quotient_coeffs.data();

    for (std::size_t i = 0; i <= dn; ++i) {
        const auto factor = rem[i] / leading;
        quo[i] = factor;

        // Вычитаем текущий член частного, умноженный на делитель
        for (std::size_t j = 0; j < divisor.n; ++j)
            rem[i + j] -= factor * divisor.c[j];
    }

    return {std::move(quotient_coeffs)};
}

Polynomial Polynomial::operator%(const Polynomial& divisor) const {
    if (divisor.isZero())
        throw std::invalid_argument("Division by zero polynomial");

    if (n < divisor.n) // Когда степень делителя больше делимого
        return {*this};

    const auto leading = divisor.c[0];
    const auto dn = n - divisor.n;

    auto remainder_coeffs = coeffs;
    const auto rem = remainder_coeffs.data();

    for (std::size_t i = 0; i <= dn; ++i) {
        const auto factor = rem[i] / leading;
        // Вычитаем текущий член частного, умноженный на делитель
        for (std::size_t j = 0; j < divisor.n; ++j)
            rem[i + j] -= factor * divisor.c[j];
    }

    // Формируем остаток (последние m коэффициентов)
    std::vector <Type> rem_coeffs;
    if (divisor.n > 1)
        rem_coeffs = std::vector<Type>(
                remainder_coeffs.begin() + int(dn + 1),
                remainder_coeffs.end()
        );

    return {std::move(rem_coeffs)};
}
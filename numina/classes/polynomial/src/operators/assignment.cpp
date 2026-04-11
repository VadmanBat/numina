//
// Created by Vadim on 18.06.2025.
//

#include "../../include/numina/polynomial.h"

Polynomial::Polynomial(const Type& constant) :
        coeffs(1, constant),
        n(1),
        c(coeffs.data())
{

}

Polynomial::Polynomial(const std::vector <Type>& coefficients) :
        coeffs(coefficients.empty() ? std::vector <Type>(1, 0) : coefficients),
        n(coeffs.size()),
        c(coeffs.data())
{
    remove_leading_zeros();
}

Polynomial::Polynomial(std::vector <Type>&& coefficients) noexcept :
        coeffs(coefficients.empty() ? std::vector <Type>(1, 0) : std::move(coefficients)),
        n(coeffs.size()),
        c(coeffs.data())
{
    remove_leading_zeros();
}

Polynomial::Polynomial(const Polynomial& other) :
        coeffs(other.coeffs),
        n(coeffs.size()),
        c(coeffs.data())
{

}

Polynomial::Polynomial(Polynomial&& other) noexcept :
        coeffs(std::move(other.coeffs)),
        n(coeffs.size()),
        c(coeffs.data())
{

}

Polynomial& Polynomial::operator=(const std::vector <Type>& coefficients) {
    coeffs = coefficients.empty() ? std::vector <Type>(1, 0) : coefficients;
    n = coeffs.size();
    c = coeffs.data();
    remove_leading_zeros();
    return *this;
}

Polynomial& Polynomial::operator=(std::vector <Type>&& coefficients) noexcept {
    coeffs = coefficients.empty() ? std::vector <Type>(1, 0) : std::move(coefficients);
    n = coeffs.size();
    c = coeffs.data();
    remove_leading_zeros();
    return *this;
}

Polynomial& Polynomial::operator=(const Polynomial& other) {
    if (&other == this)
        return *this;

    coeffs = other.coeffs;
    n = coeffs.size();
    c = coeffs.data();
    return *this;
}

Polynomial& Polynomial::operator=(Polynomial&& other) noexcept {
    coeffs = std::move(other.coeffs);
    n = coeffs.size();
    c = coeffs.data();
    return *this;
}

Polynomial& Polynomial::operator=(const Type& constant) {
    coeffs.assign(1, constant);
    n = 1;
    c = coeffs.data();
    return *this;
}
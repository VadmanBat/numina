#include "numina/poly-solver.h"
// Created by Vadim on 14.07.2025.
namespace numina {
void PolySolver::compute_derivative() const {
    const auto n = degree - 1;
    int power = static_cast<int>(degree);
    for (std::size_t i = 0; i < n; ++i) {
        c_d1[i] = power * c[i];
        c_d2[i] = --power * c_d1[i];
    }
    c_d1[n] = power * c[n];
}

void PolySolver::deflate(const Type& root, const std::size_t m) {
    if ((degree -= m) == 0)
        return;

    for (int k = 0; k < m; ++k)
        for (std::size_t i = 1; i <= degree; ++i)
            c[i] += root * c[i - 1];

    compute_derivative();
}

void PolySolver::deflate_conj(const Complex& root, const std::size_t m) {
    if ((degree -= 2 * m) == 0)
        return;

    const auto real = root.real();
    const auto imag = root.imag();
    const auto c1 = 2 * real;
    const auto c2 = real * real + imag * imag;

    for (std::size_t k = 0; k < m; ++k) {
        c[1] += c1 * c[0];
        for (std::size_t i = 2; i <= degree; ++i)
            c[i] += c1 * c[i - 1] - c2 * c[i - 2];
    }

    compute_derivative();
}
}

/* backward deflation:
void PolySolver::deflate(const Type& root, int m) {
    const auto old_degree = degree;
    if ((degree -= m) == 0)
        return;

    if (std::abs(root) > 1.0) {
        std::reverse(c, c + old_degree + 1);
        for (int k = 0; k < m; ++k)
            for (std::size_t i = 1; i <= degree; ++i)
                c[i] += c[i - 1] / root;
        std::reverse(c, c + degree + 1);
    }
    else {
        for (int k = 0; k < m; ++k)
            for (std::size_t i = 1; i <= degree; ++i)
                c[i] += root * c[i - 1];
    }

    compute_derivative();
}*/

/* compensated summation (Kahan algorithm):
void PolySolver::deflate(const Type& root, int m) {
    if ((degree -= m) == 0)
        return;

    Type err[degree + 1];
    std::fill(err, err + degree + 1, 0);

    for (int k = 0; k < m; ++k)
        for (std::size_t i = 1; i <= degree; ++i) {
            Type p = root * c[i - 1];
            Type y = p - err[i];
            Type t = c[i] + y;
            err[i] = (t - c[i]) - y;
            c[i] = t;
        }

    compute_derivative();
}

void PolySolver::deflate_conj(const Complex& root, int m) {
    if ((degree -= 2 * m) == 0)
        return;

    const auto real = root.real();
    const auto imag = root.imag();
    const auto c1 = 2 * real;
    const auto c2 = real * real + imag * imag;

    Complex err[degree + 1];
    std::fill(err, err + degree + 1, 0);

    for (int k = 0; k < m; ++k) {
        Complex p1 = c1 * c[0];
        Complex y1 = p1 - err[1];
        Complex t1 = c[1] + y1;
        err[1] = (t1 - c[1]) - y1;
        c[1] = t1;

        for (std::size_t i = 2; i <= degree; ++i) {
            Complex p = c1 * c[i - 1] - c2 * c[i - 2];
            Complex y = p - err[i];
            Complex t = c[i] + y;
            err[i] = (t - c[i]) - y;
            c[i] = t;
        }
    }

    compute_derivative();
}
*/

/* Биномиальная дефляция (более неточная)
    auto line = binom_table[m - 1];
    for (std::size_t i = degree; i >= 2; --i) {
        Type sum = 0;
        for (std::size_t j = 0; j < i; ++j)
            (sum += line[i - j] * c[j]) *= root;
        c[i] += sum;
    }
    c[1] += m * root * c[0];
*/
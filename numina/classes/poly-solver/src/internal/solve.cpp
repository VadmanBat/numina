//
// Created by Vadim on 21.08.2025.
//

#include "../../include/numina/poly-solver.h"

std::vector <PolySolver::Complex> PolySolver::solve() {
    std::vector <Complex> roots;
    roots.reserve(degree);

    Type error_previous, error_current;
    Complex x, fx, d1fx, d2fx, b, discriminant, a1, a2, gx, root;

    while (degree > 0) {
        x = Complex(-1, -1);
        for (int i = 0; i < 20; ++i) {
            fx = f(x);
            if (std::abs(fx) < e1 * std::abs(x))
                break;

            d1fx = d1(x);
            d2fx = d2(x);
            b = static_cast<Type>(degree - 1) * d1fx;
            discriminant = std::sqrt(b * b - static_cast<Type>(degree * (degree - 1)) * fx * d2fx);
            a1 = d1fx - discriminant;
            a2 = d1fx + discriminant;
            x -= static_cast<Type>(degree) * fx / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        }

        int m = 1;
        if (std::abs(fx) < e2 * std::abs(x)) {
            m = compute_multiplicity(x);
            if (m > 2)
                for (std::size_t i = dfx.size() - 1; i < m; ++i)
                    dfx.emplace_back(dfx[i].derivative());
        }

        const auto& g = dfx[m - 1];
        const auto& dg = dfx[m];

        root = x;
        gx = g(x);
        error_previous = std::abs(gx);
        for (int i = 0; i < 20; ++i) {
            x -= gx / dg(x);
            gx = g(x);
            error_current = std::abs(gx);
            if (error_current >= error_previous) {
                if (error_current == error_previous)
                    root = x;
                break;
            }
            error_previous = error_current;
            root = x;
        }

        if (std::abs(root.imag()) > e1 * std::abs(root)) {
            for (int i = 0; i < m; ++i) {
                roots.emplace_back(root);
                roots.emplace_back(std::conj(root));
            }
            deflate_conj(root, m);
        }
        else {
            const auto r = root.real();
            roots.resize(roots.size() + m, r);
            deflate(r, m);
        }
    }

    clear_state();
    return roots;
}

PolySolver::Roots PolySolver::solve_with_multiplicities() {
    Roots roots;

    Type error_previous, error_current;
    Complex x, fx, d1fx, d2fx, b, discriminant, a1, a2, gx, root;

    while (degree > 0) {
        x = Complex(-1, -1);
        for (int i = 0; i < 20; ++i) {
            fx = f(x);
            if (std::abs(fx) < e1 * std::abs(x))
                break;

            d1fx = d1(x);
            d2fx = d2(x);
            b = static_cast<Type>(degree - 1) * d1fx;
            discriminant = std::sqrt(b * b - static_cast<Type>(degree * (degree - 1)) * fx * d2fx);
            a1 = d1fx - discriminant;
            a2 = d1fx + discriminant;
            x -= static_cast<Type>(degree) * fx / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        }

        int m = 1;
        if (std::abs(fx) < e2 * std::abs(x)) {
            m = compute_multiplicity(x);
            if (m > 2)
                for (std::size_t i = dfx.size() - 1; i < m; ++i)
                    dfx.emplace_back(dfx[i].derivative());
        }

        const auto& g = dfx[m - 1];
        const auto& dg = dfx[m];

        root = x;
        gx = g(x);
        error_previous = std::abs(gx);
        for (int i = 0; i < 20; ++i) {
            x -= gx / dg(x);
            gx = g(x);
            error_current = std::abs(gx);
            if (error_current >= error_previous) {
                if (error_current == error_previous)
                    root = x;
                break;
            }
            error_previous = error_current;
            root = x;
        }

        if (std::abs(root.imag()) > e1 * std::abs(root)) {
            roots.second.emplace_back(root, m);
            deflate_conj(root, m);
        }
        else {
            const auto r = root.real();
            roots.first.emplace_back(r, m);
            deflate(r, m);
        }
    }

    clear_state();
    return roots;
}
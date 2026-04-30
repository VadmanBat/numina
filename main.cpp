#include "numina/polynomial/solve.hpp"

#include <iostream>
#include <chrono>
#include <bits/stdc++.h>

#include "numina/laplace-expression.hpp"

#include "numina/classes/polynomial/include/numina/polynomial.h"
#include "numina/classes/poly-solver/include/numina/poly-solver.h"

int main(int argc, char* argv[]) {
    const Polynomial p1({1, -2});
    const Polynomial p2({1, -6});
    const Polynomial p3({1, -10});
    const Polynomial p4({1, -20});

    //Polynomial mul(p1 * p2 * p2 * p3 * p3 * p3 * p4 * p4 * p4 * p4);
    Polynomial mul(p1 * p2 * p3 * p4 * p1 * p1 * p2 * p2 * p2 * p3 * p2 * 10000);

    /*
    std::set<int> m;
    for (int i = -1e7; i <= 1e7; ++i) {
        auto x = 2 + i / 1e8; // 2, 1e6, 1e8
        if (auto v = mul.computeMultiplicity(x); v == 3)
            m.insert(v);
        else {
            std::cout << "m(" << x << ") = " << v << '\n';
            for (int i = 1; i <= 5; ++i) {
                std::cout << std::fixed << std::setprecision(20);
                std::cout << "m = " << i << ", f(x, m) = " << mul.multiplicity(x, i) << '\n';
            }
        }
        //mul.computeMultiplicity(x);
    }
    std::cout << m.size() << '\n';
    for (auto mmm : m)
        std::cout << mmm << ' ';
    std::cout << '\n';
    return 0;*/
    //return 0;
    //Polynomial mul(p1 * p2 * p3);

    //Polynomial mul(p1 * p2 * p3 * p4);
    //Polynomial mul(p1 * p2 * p3 * p4);
    //Polynomial mul(Polynomial({1, -5}) * Polynomial({1, -2, 1.01}));

    Polynomial sum(mul + p1);

    /*std::cout << mul << '\n';
    std::cout << sum << '\n';
    std::cout << gcd(mul, mul.derivative()) << '\n';*/

    std::cout << "poly-solver:\n";
    std::cout << std::fixed << std::setprecision(20);
    for (auto root : mul.computeRoots())
        std::cout << root << '\n';
    std::cout << '\n';

    std::cout << "poly-solver-2:\n";
    std::cout << std::fixed << std::setprecision(20);
    numina::PolySolver x;
    auto roots = x.multisolve(mul);
    for (auto [root, m] : roots.first)
        std::cout << root << ' ' << m << '\n';
    for (auto [root, m] : roots.second)
        std::cout << root << ' ' << m << '\n';
    std::cout << '\n';

    std::cout << "func-solver:\n";
    for (auto root : numina::solve_polynomial_laguerre(mul.vector()))
        std::cout << root << '\n';
    std::cout << '\n';

    {
        auto start = std::chrono::high_resolution_clock::now();

        std::vector coeffs = mul.vector();
        for (int i = 0; i < 1e5; ++i)
            auto roots = numina::solve_polynomial_laguerre(coeffs);

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "time: " << (end - start).count() / 1e9 << '\n';
        std::cout << "sum = " << sum << '\n';
    }
    {
        auto start = std::chrono::high_resolution_clock::now();

        numina::PolySolver x;
        for (int i = 0; i < 1e5; ++i) {
            //std::cout << i << '\n';
            auto roots = x.multisolve(mul, numina::PolySolver::Method::Explicit);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "time: " << (end - start).count() / 1e9 << '\n';
        std::cout << "sum = " << sum << '\n';
    }

    {
        auto start = std::chrono::high_resolution_clock::now();

        numina::PolySolver x;
        for (int i = 0; i < 1e5; ++i) {
            //std::cout << i << '\n';
            auto roots = x.solve(mul, numina::PolySolver::Method::Implicit);
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "time: " << (end - start).count() / 1e9 << '\n';
        std::cout << "sum = " << sum << '\n';
    }

    /*TranFunc tf;

    test(std::bind(&TermExpression<double>::operator(), &tf.impulseResp, std::placeholders::_1));
    //test(std::bind(&TranFunc::impulseResponse, &tf, std::placeholders::_1));
    //test(std::bind(&TranFunc::fullImpulseResponse, &tf, std::placeholders::_1));

    return 0;*/
}
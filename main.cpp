#include "numina/terms/term-expression.hpp"
#include "numina/polynomial/solve.hpp"

#include <iostream>
#include <chrono>
#include <functional>
#include <bits/stdc++.h>

class TranFunc {
private:
    std::vector<std::complex<double>> roots = {
        {-0.05, -0.5}, {-0.05, 0.5}, {-0.05, -0.5}, {-0.05, 0.5},
        {-0.05, -0.5}, {-0.05, 0.5}, {-0.05, -0.5}, {-0.05, 0.5},
        {-0.05, -0.5}, {-0.05, 0.5}, {-0.05, -0.5}, {-0.05, 0.5},
        {-0.05, -0.5}, {-0.05, 0.5}, {-0.05, -0.5}, {-0.05, 0.5}
    };

    //std::vector <Term <std::complex<double>>*> terms;
    std::vector<Term<double>*> terms;

    std::vector<int> powers = {0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2};
    //std::vector <double> powers = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    //std::vector <double> coeffs = {5, 4, 8, 1, 5, 4, 8, 1, 5, 4, 8, 1, 5, 4, 8, 1};
    std::vector<std::complex<double>> coeffs = {5, 4, 8, 1, 5, 4, 8, 1, 5, 4, 8, 1, 5, 4, 8, 1};

public:
    TranFunc() : impulseResp(roots, coeffs, powers) {
        /*const std::size_t n = roots.size();
        for (std::size_t i = 0; i < n; i += 2) {
            terms.push_back(createTerm(coeffs[i], roots[i], (int)powers[i]));
        }

        for (auto t : terms)
            std::cout << t->string() << '\n';*/
    }

    TermExpression<double> impulseResp;

    ~TranFunc() {
        for (auto term : terms)
            delete term;
    }

    inline double impulseResponse(double time) {
        std::complex<double> answer = 0;
        const std::size_t n         = roots.size();
        for (std::size_t i = 0; i < n; ++i)
            answer += coeffs[i] * std::exp(roots[i] * time);

        return answer.real();
    }

    inline double fullImpulseResponse(double time) {
        std::complex<double> answer = 0;
        const std::size_t n         = roots.size();
        for (std::size_t i = 0; i < n; ++i)
            answer += coeffs[i] * std::exp(roots[i]) * std::pow(time, powers[i]);

        return answer.real();
    }

    inline double altImpulseResponse(double time) {
        //std::complex <double> answer = 0;
        double answer = 0;
        //Term<std::complex<double>>::time = time;
        Term<double>::time  = time;
        const std::size_t n = terms.size();
        for (std::size_t i = 0; i < n; ++i)
            answer += terms[i]->value();

        return answer;
    }
};

void test(const std::function<double(double)>& h) {
    auto start = std::chrono::high_resolution_clock::now();

    double sum = 0;
    for (int i = 0; i < 1e7; ++i)
        sum += h(i);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "time: " << (end - start).count() / 1e9 << '\n';
    std::cout << "sum = " << sum << '\n';
}

#include "numina/classes/polynomial/include/numina/polynomial.h"
#include "numina/classes/poly-solver/include/numina/poly-solver.h"

int main(int argc, char* argv[]) {
    const Polynomial p1({1, -2});
    const Polynomial p2({1, -6});
    const Polynomial p3({1, -10});
    const Polynomial p4({1, -20});

    //Polynomial mul(p1 * p2 * p2 * p3 * p3 * p3 * p4 * p4 * p4 * p4);
    Polynomial mul(p1 * p2 * p3 * p4 * p1 * p1 * p2 * p2 * p2 * p3 * p2);

    std::set<int> m;
    for (int i = -1e1; i <= 1e1; ++i) {
        auto x = 6 + i / 1e7; // 2, 1e6, 1e8
        if (auto v = mul.__computeMultiplicity(x); v == 1+1)
            m.insert(v);
        else {
            std::cout << "m(" << x << ") = " << v << '\n';
            auto mmm = mul.__computeMultiplicity(x, true);

            //std::cout << "x = " << x << ", m = " << mul.__computeMultiplicity(x) << '\n';
            for (int i = 1; i <= 5; ++i) {
                std::cout << std::fixed << std::setprecision(20);
                std::cout << "m = " << i << ", f(x, m) = " << mul.multiplicity(x, i) << '\n';
                std::cout << "m = " << i << ", f(x, m) = " << mul.multiplicity((long double)(x), i) << '\n';
                std::cout << "m = " << i << ", f(x, m) = " << mul.lag(x, i) << '\n';
                std::cout << "m = " << i << ", f(x, m) = " << mul.lag((long double)(x), i) << '\n';
            }
        }
        //mul.computeMultiplicity(x);
    }
    std::cout << m.size() << '\n';
    for (auto mmm : m)
        std::cout << mmm << ' ';
    std::cout << '\n';
    return 0;
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
    auto roots = x.solveWithImplicitDeflation(mul);
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
            auto roots = x.solveWithImplicitDeflation(mul);
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
            auto roots = x.solve(mul);
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
#pragma once
// Created by Vadim on 14.04.2026.
#include <map>

namespace numina {
class QuadraticDerivPoly {
public:
    using ll    = long long;
    using Terms = std::map<std::pair<int, int>, ll>;

private:
    Terms terms; // i <= j : coeff при p^{(i)} * p^{(j)}

public:
    QuadraticDerivPoly() = default;
    QuadraticDerivPoly(std::initializer_list<std::tuple<int, int, ll>> monomials);

    void addTerm(int i, int j, ll coeff);
    const Terms& getTerms() const noexcept;

    void differentiate();       // d/dt (p^{(a)} p^{(b)}) = p^{(a+1)} p^{(b)} + p^{(a)} p^{(b+1)}
    void doubleDifferentiate(); // d²/dt² (p^{(a)} p^{(b)}) = p^{(a+2)} p^{(b)} + 2 * p^(a+1)p^(b+1) + p^{(a)} p^{(b+2)}
};
}
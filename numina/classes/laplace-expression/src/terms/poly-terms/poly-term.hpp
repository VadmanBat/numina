#pragma once
// Created by Vadim on 07.01.2025.
#include "../time-terms/time-term.hpp"
#include <cmath>

namespace numina {
template <typename Type>
class PolyTerm : public Term<Type> {
    const Type coefficient;
    const int power;

public:
    explicit PolyTerm(Type a, int n) :
        coefficient(a),
        power(n) {
    }

    PolyTerm(const PolyTerm& other) :
        coefficient(other.coefficient),
        power(other.power) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<PolyTerm>(*this);
    }

    Type value(Type t) const override {
        return coefficient * std::pow(t, power);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        if (power == 2)
            return {std::make_unique<TimeTerm>(2 * coefficient)};
        return {std::make_unique<PolyTerm>(power * coefficient, power - 1)};
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] std::string string() const override {
        return (std::stringstream() << coefficient << " × t<sup>" << power << "</sup>").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × t<sup>" << power << "</sup>").str();
    }
};
}
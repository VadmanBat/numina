#pragma once
// Created by Vadim on 07.01.2025.
#include "../time-terms/exp-time-term.hpp"

namespace numina {
template <typename Type>
class ExpPolyTerm : public Term<Type> {
    const Type coefficient, root;
    const int power;

public:
    explicit ExpPolyTerm(Type c, Type r, int n) :
        coefficient(c),
        root(r),
        power(n) {
    }

    ExpPolyTerm(const ExpPolyTerm& other) :
        coefficient(other.coefficient),
        root(other.root),
        power(other.power) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<ExpPolyTerm>(*this);
    }

    Type value(Type t) const override {
        return coefficient * std::exp(root * t) * std::pow(t, power);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        std::vector<std::unique_ptr<Term<Type>>> result;
        result.reserve(2);
        if (power == 2)
            result.push_back(std::make_unique<ExpTimeTerm<Type>>(2 * coefficient, root));
        else
            result.push_back(std::make_unique<ExpPolyTerm>(power * coefficient, root, power - 1));
        result.push_back(std::make_unique<ExpPolyTerm>(root * coefficient, root, power));
        return result;
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup> × t<sup>" << power << "</sup>")
            .str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup> × t<sup>" << power <<
                "</sup>").str();
    }
};
}
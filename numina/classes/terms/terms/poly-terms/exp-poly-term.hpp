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

    inline Term<Type>* clone() const override {
        return new ExpPolyTerm(*this);
    }

    inline Type value() const override {
        return coefficient * std::exp(root * Term<Type>::time) * std::pow(Term<Type>::time, power);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        if (power == 2)
            return {
                new ExpTimeTerm(2 * coefficient, root),
                new ExpPolyTerm(root * coefficient, root, 2)
            };
        return {
            new ExpPolyTerm(power * coefficient, root, power - 1),
            new ExpPolyTerm(root * coefficient, root, power)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup> × t<sup>" << power << "</sup>")
            .str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup> × t<sup>" << power <<
                "</sup>").str();
    }
};
}
#pragma once
// Created by Vadim on 07.01.2025.
#include "../terms/exp-term.hpp"

namespace numina {
template <typename Type>
class ExpTimeTerm : public Term<Type> {
private:
    Type coefficient, root;

public:
    explicit ExpTimeTerm(Type c, Type r) :
        coefficient(c),
        root(r) {
    }

    ExpTimeTerm(const ExpTimeTerm& other) :
        coefficient(other.coefficient),
        root(other.root) {
    }

    inline Term<Type>* clone() const override {
        return new ExpTimeTerm(*this);
    }

    inline Type value() const override {
        return coefficient * std::exp(root * Term<Type>::time) * Term<Type>::time;
    }

    inline std::vector<Term<Type>*> derivative() const override {
        return {
            new ExpTerm(coefficient * root, root),
            new ExpTerm(coefficient, root)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup> × t").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup> × t").str();
    }
};
}
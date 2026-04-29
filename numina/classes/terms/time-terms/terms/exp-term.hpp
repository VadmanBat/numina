#pragma once
// Created by Vadim on 05.01.2025.
#include "term.hpp"
#include <cmath>

namespace numina {
template <typename Type>
class ExpTerm : public Term<Type> {
private:
    const Type coefficient, root;

public:
    explicit ExpTerm(Type c, Type r) :
        coefficient(c),
        root(r) {
    }

    ExpTerm(const ExpTerm& other) :
        coefficient(other.coefficient),
        root(other.root) {
    }

    inline Term<Type>* clone() const override {
        return new ExpTerm(*this);
    }

    inline Type value() const override {
        return coefficient * std::exp(root * Term<Type>::time);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        return {new ExpTerm(coefficient * root, root)};
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup>").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup>").str();
    }
};
}
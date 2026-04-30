#pragma once
// Created by Vadim on 07.01.2025.
#include "../base-terms/exp-term.hpp"

namespace numina {
template <typename Type>
class ExpTimeTerm : public Term<Type> {
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

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<ExpTimeTerm>(*this);
    }

    Type value(Type t) const override {
        return coefficient * std::exp(root * t) * t;
    }

    std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        std::vector<std::unique_ptr<Term<Type>>> result;
        result.reserve(3);
        result.push_back(std::make_unique<ExpTerm<Type>>(coefficient * root, root));
        result.push_back(std::make_unique<ExpTerm<Type>>(coefficient, root));
        return result;
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup> × t").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup> × t").str();
    }
};
}
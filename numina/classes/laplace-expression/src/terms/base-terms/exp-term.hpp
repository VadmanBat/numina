#pragma once
// Created by Vadim on 05.01.2025.
#include "../../term.hpp"
#include <cmath>

namespace numina {
template <typename Type>
class ExpTerm : public Term<Type> {
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

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<ExpTerm>(*this);
    }

    Type value(Type t) const override {
        return coefficient * std::exp(root * t);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        std::vector<std::unique_ptr<Term<Type>>> result;
        result.reserve(1);
        result.emplace_back(std::make_unique<ExpTerm>(coefficient * root, root));
        return result;
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup>").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup>").str();
    }
};
}
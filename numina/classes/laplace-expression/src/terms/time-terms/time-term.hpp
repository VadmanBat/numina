#pragma once
// Created by Vadim on 07.01.2025.
#include "../../term.hpp"

namespace numina {
template <typename Type>
class TimeTerm : public Term<Type> {
    const Type coefficient;

public:
    explicit TimeTerm(Type a) :
        coefficient(a) {
    }

    TimeTerm(const TimeTerm& other) :
        coefficient(other.coefficient) {
    }

    Term<Type>* clone() const override {
        return new TimeTerm(*this);
    }

    Type value() const override {
        return coefficient * Term<Type>::time;
    }

    [[nodiscard]] std::vector<Term<Type>*> derivative() const override {
        return {};
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] std::string string() const override {
        return (std::stringstream() << coefficient << " × t").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × t").str();
    }

    [[nodiscard]] Type derivativeConstant() const override {
        return coefficient;
    };
};
}
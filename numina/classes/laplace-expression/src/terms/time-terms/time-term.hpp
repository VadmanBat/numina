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

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<TimeTerm>(*this);
    }

    Type value(Type t) const override {
        return coefficient * t;
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
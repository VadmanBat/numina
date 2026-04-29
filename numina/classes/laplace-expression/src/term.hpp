#pragma once
// Created by Vadim on 05.01.2025.
#include <memory>
#include <vector>
#include <sstream>
#include <iomanip>

namespace numina {
template <typename Type>
class Term {
public:
    virtual ~Term() = default;
    [[nodiscard]] virtual Type value(Type t) const = 0;
    virtual std::vector<std::unique_ptr<Term>> derivative() const { return {}; }
    virtual std::unique_ptr<Term> clone() const = 0;

    [[nodiscard]] virtual bool isPositive() const = 0;
    [[nodiscard]] virtual std::string string() const = 0;
    [[nodiscard]] virtual std::string unsignedString() const = 0;

    [[nodiscard]] virtual Type derivativeConstant() const { return 0; };

    static void setPrecision(const int n) {
        std::stringstream() << std::fixed << std::setprecision(n);
    }
};
}
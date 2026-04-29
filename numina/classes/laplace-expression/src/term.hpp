#pragma once
// Created by Vadim on 05.01.2025.
#include <sstream>
#include <iomanip>
#include <vector>

namespace numina {
template <typename Type>
class Term {
public:
    virtual ~Term() = default;
    [[nodiscard]] virtual Type value(Type t) const = 0;
    [[nodiscard]] virtual std::vector<Term*> derivative() const = 0;
    [[nodiscard]] virtual Term* clone() const = 0;

    [[nodiscard]] virtual bool isPositive() const = 0;
    [[nodiscard]] virtual std::string string() const = 0;
    [[nodiscard]] virtual std::string unsignedString() const = 0;

    [[nodiscard]] virtual Type derivativeConstant() const { return 0; };

    static void setPrecision(const int n) {
        std::stringstream() << std::fixed << std::setprecision(n);
    }
};
}
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
    [[nodiscard]] virtual inline Type value() const = 0;
    [[nodiscard]] virtual inline std::vector<Term*> derivative() const = 0;
    [[nodiscard]] virtual inline Term* clone() const = 0;

    [[nodiscard]] virtual inline bool isPositive() const = 0;
    [[nodiscard]] virtual inline std::string string() const = 0;
    [[nodiscard]] virtual inline std::string unsignedString() const = 0;

    [[nodiscard]] virtual inline Type derivativeConstant() const {
        return 0;
    };

    inline static Type time;

    static void setPrecision(int n) {
        std::stringstream() << std::fixed << std::setprecision(n);
    }
};
}
#pragma once
// Created by Vadim on 05.01.2025.
#include "../../term.hpp"
#include <complex>

namespace numina {
template <typename Type>
class CosTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    const Type amplitude, omega, phi;

public:
    explicit CosTerm(const Comp& c, const Comp& r) :
        amplitude(2 * std::abs(c)),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())) {
    }

    explicit CosTerm(Type amplitude, Type omega, Type phi) :
        amplitude(amplitude),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)) {
    }

    CosTerm(const CosTerm& other) :
        amplitude(other.amplitude),
        omega(other.omega),
        phi(other.phi) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<CosTerm>(*this);
    }

    Type value(Type t) const override {
        return amplitude * std::cos(omega * t + phi);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        std::vector<std::unique_ptr<Term<Type>>> result;
        result.reserve(1);
        result.push_back(std::make_unique<CosTerm>(-amplitude * omega, omega, phi - std::numbers::pi_v<Type> / 2));
        return result;
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t)").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ')').str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t)").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ')').str();
    }
};
}
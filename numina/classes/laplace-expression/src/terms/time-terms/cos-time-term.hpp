#pragma once
// Created by Vadim on 07.01.2025.
#include "../base-terms/cos-term.hpp"

namespace numina {
template <typename Type>
class CosTimeTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    const Type amplitude, omega, phi;

public:
    explicit CosTimeTerm(const Comp& c, const Comp& r) :
        amplitude(2 * std::abs(c)),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())) {
    }

    explicit CosTimeTerm(Type amplitude, Type omega, Type phi) :
        amplitude(amplitude),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)) {
    }

    CosTimeTerm(const CosTimeTerm& other) :
        amplitude(other.amplitude),
        omega(other.omega),
        phi(other.phi) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<CosTimeTerm>(*this);
    }

    Type value(Type t) const override {
        return amplitude * std::cos(omega * t + phi) * t;
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        return {
            std::make_unique<CosTerm<Type>>(amplitude, omega, phi),
            std::make_unique<CosTimeTerm>(-amplitude * omega, omega, phi - std::numbers::pi_v<Type> / 2)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t").str();
    }
};
}
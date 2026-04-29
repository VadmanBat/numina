#pragma once
// Created by Vadim on 07.01.2025.
#include "../base-terms/exp-cos-term.hpp"

namespace numina {
template <typename Type>
class ExpCosTimeTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    const Type amplitude, alpha, omega, phi;

public:
    explicit ExpCosTimeTerm(const Comp& c, const Comp& r) :
        amplitude(2 * std::abs(c)),
        alpha(r.real()),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())) {
    }

    explicit ExpCosTimeTerm(Type amplitude, Type alpha, Type omega, Type phi) :
        amplitude(amplitude),
        alpha(alpha),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)) {
    }

    ExpCosTimeTerm(const ExpCosTimeTerm& other) :
        amplitude(other.amplitude),
        alpha(other.alpha),
        omega(other.omega),
        phi(other.phi) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<ExpCosTimeTerm>(*this);
    }

    Type value(Type t) const override {
        return amplitude * std::exp(alpha * t) * std::cos(omega * t + phi) * t;
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        return {
            std::make_unique<ExpCosTerm<Type>>(amplitude, alpha, omega, phi),
            std::make_unique<ExpCosTimeTerm>(amplitude * alpha, alpha, omega, phi),
            std::make_unique<ExpCosTimeTerm>(-amplitude * omega, alpha, omega, phi - std::numbers::pi_v<Type> / 2)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ") × t").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ") × t").str();
    }
};
}
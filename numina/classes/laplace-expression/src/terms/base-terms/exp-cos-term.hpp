#pragma once
// Created by Vadim on 05.01.2025.
#include "../../term.hpp"
#include <complex>

namespace numina {
template <typename Type>
class ExpCosTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    const Type amplitude, alpha, omega, phi;

public:
    explicit ExpCosTerm(const Comp& c, const Comp& r) :
        amplitude(2 * std::abs(c)),
        alpha(r.real()),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())) {
    }

    explicit ExpCosTerm(Type amplitude, Type alpha, Type omega, Type phi) :
        amplitude(amplitude),
        alpha(alpha),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)) {
    }

    ExpCosTerm(const ExpCosTerm& other) :
        amplitude(other.amplitude),
        alpha(other.alpha),
        omega(other.omega),
        phi(other.phi) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return std::make_unique<ExpCosTerm>(*this);
    }

    Type value(Type t) const override {
        return amplitude * std::exp(alpha * t) * std::cos(omega * t + phi);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        std::vector<std::unique_ptr<Term<Type>>> result;
        result.reserve(2);
        result.push_back(std::make_unique<ExpCosTerm>(amplitude * alpha, alpha, omega, phi));
        result.push_back(std::make_unique<ExpCosTerm>(-omega * amplitude, alpha, omega, phi - std::numbers::pi_v<Type> / 2));
        return result;
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t)")
                .str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ')').str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t)")
                .str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ')').str();
    }
};
}
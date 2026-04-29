#pragma once
// Created by Vadim on 07.01.2025.
#include "../time-terms/exp-cos-time-term.hpp"

namespace numina {
template <typename Type>
class ExpCosPolyTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    Type amplitude, alpha, omega, phi;
    int power;

public:
    explicit ExpCosPolyTerm(const Comp& c, const Comp& r, int n) :
        amplitude(2 * std::abs(c)),
        alpha(r.real()),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())),
        power(n) {
    }

    explicit ExpCosPolyTerm(Type amplitude, Type alpha, Type omega, Type phi, int power) :
        amplitude(amplitude),
        alpha(alpha),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)),
        power(power) {
    }

    ExpCosPolyTerm(const ExpCosPolyTerm& other) :
        amplitude(other.amplitude),
        alpha(other.alpha),
        omega(other.omega),
        phi(other.phi),
        power(other.power) {
    }

    std::unique_ptr<Term<Type>> clone() const override {
        return new ExpCosPolyTerm(*this);
    }

    Type value(Type t) const override {
        return amplitude * std::exp(alpha * t) * std::cos(omega * t + phi) * std::pow(t, power);
    }

    [[nodiscard]] std::vector<std::unique_ptr<Term<Type>>> derivative() const override {
        if (power == 2)
            return {
                std::make_unique<ExpCosTimeTerm>(2 * amplitude, alpha, omega, phi),
                std::make_unique<ExpCosPolyTerm>(alpha * amplitude, alpha, omega, phi, 2),
                std::make_unique<ExpCosPolyTerm>(-amplitude * omega, alpha, omega, phi - std::numbers::pi_v<Type> / 2, 2)
            };
        return {
            std::make_unique<ExpCosPolyTerm>(power * amplitude, alpha, omega, phi, power - 1),
            std::make_unique<ExpCosPolyTerm>(alpha * amplitude, alpha, omega, phi, power),
            std::make_unique<ExpCosPolyTerm>(-amplitude * omega, alpha, omega, phi - std::numbers::pi_v<Type> / 2, power)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t<sup>" << power << "</sup>").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << " × t<sup>" << power << "</sup>").str();
    }

    [[nodiscard]] std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t<sup>" << power << "</sup>").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << " × t<sup>" << power << "</sup>").str();
    }
};
}
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

    inline Term<Type>* clone() const override {
        return new ExpCosTerm(*this);
    }

    inline Type value() const override {
        return amplitude * std::exp(alpha * Term<Type>::time) * std::cos(omega * Term<Type>::time + phi);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        return {
            new ExpCosTerm(amplitude * alpha, alpha, omega, phi),
            new ExpCosTerm(-omega * amplitude, alpha, omega, phi - std::numbers::pi_v<Type> / 2)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t)")
                .str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ')').str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t)")
                .str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << ')').str();
    }
};
}
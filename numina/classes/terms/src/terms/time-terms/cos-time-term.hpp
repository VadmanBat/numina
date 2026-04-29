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

    inline Term<Type>* clone() const override {
        return new CosTimeTerm(*this);
    }

    inline Type value() const override {
        return amplitude * std::cos(omega * Term<Type>::time + phi) * Term<Type>::time;
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        return {
            new CosTerm(amplitude, omega, phi),
            new CosTimeTerm(-amplitude * omega, omega, phi - std::numbers::pi_v<Type> / 2)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t").str();
    }
};
}
#pragma once
// Created by Vadim on 07.01.2025.
#include "../time-terms/cos-time-term.hpp"

namespace numina {
template <typename Type>
class CosPolyTerm : public Term<Type> {
    using Comp = std::complex<Type>;
    const Type amplitude, omega, phi;
    const int power;

public:
    explicit CosPolyTerm(const Comp& c, const Comp& r, int n) :
        amplitude(2 * std::abs(c)),
        omega(r.imag()),
        phi(std::atan2(c.imag(), c.real())),
        power(n) {
    }

    explicit CosPolyTerm(Type amplitude, Type omega, Type phi, int power) :
        amplitude(amplitude),
        omega(omega),
        phi(std::remainder(phi, 2 * std::numbers::pi_v<Type>)),
        power(power) {
    }

    CosPolyTerm(const CosPolyTerm& other) :
        amplitude(other.amplitude),
        omega(other.omega),
        phi(other.phi),
        power(other.power) {
    }

    inline Term<Type>* clone() const override {
        return new CosPolyTerm(*this);
    }

    inline Type value() const override {
        return amplitude * std::cos(omega * Term<Type>::time + phi) * std::pow(Term<Type>::time, power);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        if (power == 2)
            return {
                new CosTimeTerm(2 * amplitude, omega, phi),
                new CosPolyTerm(-amplitude * omega, omega, phi - std::numbers::pi_v<Type> / 2, 2)
            };
        return {
            new CosPolyTerm(power * amplitude, omega, phi, power - 1),
            new CosPolyTerm(-amplitude * omega, omega, phi - std::numbers::pi_v<Type> / 2, power)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t<sup>" << power << "</sup>").
                str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t<sup>" << power << "</sup>").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × cos(" << omega << " × t) × t<sup>" << power << "</sup>").
                str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × cos(" << omega << " × t " << phi_sign << ' ' << std::abs(phi) <<
                ") × t<sup>" << power << "</sup>").str();
    }
};
}
#pragma once
// Created by Vadim on 07.01.2025.
#include "time-terms.hpp"

namespace numina {
template <typename Type>
class PolyTerm : public Term<Type> {
private:
    const Type coefficient;
    const int power;

public:
    explicit PolyTerm(Type a, int n) :
        coefficient(a),
        power(n) {
    }

    PolyTerm(const PolyTerm& other) :
        coefficient(other.coefficient),
        power(other.power) {
    }

    inline Term<Type>* clone() const override {
        return new PolyTerm(*this);
    }

    inline Type value() const override {
        return coefficient * std::pow(Term<Type>::time, power);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        if (power == 2)
            return {new TimeTerm(2 * coefficient)};
        return {new PolyTerm(power * coefficient, power - 1)};
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        return (std::stringstream() << coefficient << " × t<sup>" << power << "</sup>").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × t<sup>" << power << "</sup>").str();
    }
};

template <typename Type>
class ExpPolyTerm : public Term<Type> {
private:
    const Type coefficient, root;
    const int power;

public:
    explicit ExpPolyTerm(Type c, Type r, int n) :
        coefficient(c),
        root(r),
        power(n) {
    }

    ExpPolyTerm(const ExpPolyTerm& other) :
        coefficient(other.coefficient),
        root(other.root),
        power(other.power) {
    }

    inline Term<Type>* clone() const override {
        return new ExpPolyTerm(*this);
    }

    inline Type value() const override {
        return coefficient * std::exp(root * Term<Type>::time) * std::pow(Term<Type>::time, power);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        if (power == 2)
            return {
                new ExpTimeTerm(2 * coefficient, root),
                new ExpPolyTerm(root * coefficient, root, 2)
            };
        return {
            new ExpPolyTerm(power * coefficient, root, power - 1),
            new ExpPolyTerm(root * coefficient, root, power)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return coefficient > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        return (std::stringstream() << coefficient << " × e<sup>" << root << " × t</sup> × t<sup>" << power << "</sup>")
            .str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        return (std::stringstream() << std::abs(coefficient) << " × e<sup>" << root << " × t</sup> × t<sup>" << power <<
                "</sup>").str();
    }
};

template <typename Type>
class CosPolyTerm : public Term<Type> {
private:
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

template <typename Type>
class ExpCosPolyTerm : public Term<Type> {
private:
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

    inline Term<Type>* clone() const override {
        return new ExpCosPolyTerm(*this);
    }

    inline Type value() const override {
        return amplitude * std::exp(alpha * Term<Type>::time) * std::cos(omega * Term<Type>::time + phi) * std::pow(
                   Term<Type>::time, power);
    }

    [[nodiscard]] inline std::vector<Term<Type>*> derivative() const override {
        if (power == 2)
            return {
                new ExpCosTimeTerm(2 * amplitude, alpha, omega, phi),
                new ExpCosPolyTerm(alpha * amplitude, alpha, omega, phi, 2),
                new ExpCosPolyTerm(-amplitude * omega, alpha, omega, phi - std::numbers::pi_v<Type> / 2, 2)
            };
        return {
            new ExpCosPolyTerm(power * amplitude, alpha, omega, phi, power - 1),
            new ExpCosPolyTerm(alpha * amplitude, alpha, omega, phi, power),
            new ExpCosPolyTerm(-amplitude * omega, alpha, omega, phi - std::numbers::pi_v<Type> / 2, power)
        };
    }

    [[nodiscard]] bool isPositive() const override {
        return amplitude > 0;
    }

    [[nodiscard]] inline std::string string() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t<sup>" << power << "</sup>").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << " × t<sup>" << power << "</sup>").str();
    }

    [[nodiscard]] inline std::string unsignedString() const override {
        if (phi == 0)
            return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega <<
                    " × t) × t<sup>" << power << "</sup>").str();
        const char phi_sign = phi > 0 ? '+' : '-';
        return (std::stringstream() << amplitude << " × e<sup>" << alpha << " × t</sup> × cos(" << omega << " × t " <<
                phi_sign << ' ' << std::abs(phi) << " × t<sup>" << power << "</sup>").str();
    }
};
}
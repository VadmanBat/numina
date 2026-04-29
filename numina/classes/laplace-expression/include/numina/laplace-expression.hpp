#pragma once
// Created by Vadim on 07.01.2025.
#include "../../src/terms/poly-terms/poly-term.hpp"
#include "../../src/terms/poly-terms/exp-poly-term.hpp"
#include "../../src/terms/poly-terms/cos-poly-term.hpp"
#include "../../src/terms/poly-terms/exp-cos-poly-term.hpp"

namespace numina {
template <typename Type>
class LaplaceExpression {
    using Comp    = std::complex<Type>;
    using VecComp = std::vector<Comp>;

    Type init_value;
    std::vector<Term<Type>*> terms;

public:
    explicit LaplaceExpression() : init_value(0) {
    }

    explicit LaplaceExpression(const VecComp& roots, const VecComp& coeffs, const std::vector<int>& powers) :
        init_value(0) {
        static const Type epsilon = std::sqrt(std::numeric_limits<Type>::epsilon());
        //const auto n = std::min({roots.size(), coeffs.size(), powers.size()});
        const auto n = std::min(std::min(roots.size(), coeffs.size()), powers.size());
        for (std::size_t i = 0; i < n; ++i) {
            if (std::abs(coeffs[i]) < epsilon)
                continue;
            emplace_back(coeffs[i], roots[i], powers[i]);
            if (i + 1 < n && std::abs(roots[i] - std::conj(roots[i + 1])) < epsilon && std::abs(
                    coeffs[i] - std::conj(coeffs[i + 1])) < epsilon)
                ++i;
        }
    }

    explicit LaplaceExpression(const std::vector<Term<Type>*>& terms, Type init_value = 0) :
        init_value(init_value),
        terms(terms) {
    }

    explicit LaplaceExpression(std::vector<Term<Type>*>&& terms, Type init_value = 0) :
        init_value(init_value),
        terms(std::move(terms)) {
    }

    ~LaplaceExpression() {
        for (const auto term : terms)
            delete term;
    }

    void emplace_back(const Comp& c, const Comp& r, int n) {
        if (std::abs(r) == 0)
            switch (n) {
                case 0:
                    init_value += c.real();
                    return;
                case 1:
                    terms.push_back(new TimeTerm(c.real()));
                    return;
                default:
                    terms.push_back(new PolyTerm(c.real(), n));
                    return;
            }

        if (r.imag() == 0)
            switch (n) {
                case 0:
                    terms.push_back(new ExpTerm(c.real(), r.real()));
                    return;
                case 1:
                    terms.push_back(new ExpTimeTerm(c.real(), r.real()));
                    return;
                default:
                    terms.push_back(new ExpPolyTerm(c.real(), r.real(), n));
                    return;
            }

        if (r.real() == 0)
            switch (n) {
                case 0:
                    terms.push_back(new CosTerm(c, r));
                    return;
                case 1:
                    terms.push_back(new CosTimeTerm(c, r));
                    return;
                default:
                    terms.push_back(new CosPolyTerm(c, r, n));
                    return;
            }

        switch (n) {
            case 0:
                terms.push_back(new ExpCosTerm(c, r));
                return;
            case 1:
                terms.push_back(new ExpCosTimeTerm(c, r));
                return;
            default:
                terms.push_back(new ExpCosPolyTerm(c, r, n));
                return;
        }
    }

    void push_back(std::vector<Term<Type>*> new_terms) {
        terms.insert(terms.end(), new_terms.begin(), new_terms.end());
    }

    Type operator()(double time) const {
        Term<Type>::time = time;

        Type result = init_value;
        for (const auto term : terms)
            result += term->value();

        return result;
    }

    LaplaceExpression derivative() const {
        Type derivative_init_value = 0;
        std::vector<Term<Type>*> derivative_terms;
        for (const auto term : terms) {
            derivative_init_value += term->derivativeConstant();
            auto derivative       = term->derivative();
            derivative_terms.insert(
                derivative_terms.end(),
                derivative.begin(), derivative.end()
                );
        }
        return LaplaceExpression(derivative_terms, derivative_init_value);
    }

    [[nodiscard]] std::string string() const {
        if (init_value != 0) {
            std::string result = (std::stringstream() << init_value).str();
            for (const auto term : terms)
                result += (term->isPositive() ? " + " : " - ") + term->unsignedString();
            return result;
        }
        const auto n = terms.size();
        if (n == 0)
            return {};
        std::string result = terms[0]->string();
        for (std::size_t i = 1; i < n; ++i)
            result += (terms[i]->isPositive() ? " + " : " - ") + terms[i]->unsignedString();
        return result;
    }

    LaplaceExpression& operator=(const LaplaceExpression& other) {
        for (const auto term : terms)
            delete term;
        terms.clear();
        init_value = other.init_value;
        terms.reserve(other.terms.size());
        for (auto term : other.terms)
            terms.push_back(term->clone());
        return *this;
    }

    LaplaceExpression& operator=(LaplaceExpression&& other) noexcept {
        for (const auto term : terms)
            delete term;
        terms.clear();
        init_value = other.init_value;
        terms      = std::move(other.terms);
        return *this;
    }
};
}
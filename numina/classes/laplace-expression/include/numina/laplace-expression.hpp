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

    Type init_value = 0;
    std::vector<std::unique_ptr<Term<Type>>> terms;

public:
    explicit LaplaceExpression() = default;

    explicit LaplaceExpression(const VecComp& roots, const VecComp& coeffs, const std::vector<int>& powers) {
        static const Type epsilon = std::sqrt(std::numeric_limits<Type>::epsilon());
        const auto n = std::min({roots.size(), coeffs.size(), powers.size()});
        for (std::size_t i = 0; i < n; ++i) {
            if (std::abs(coeffs[i]) < epsilon)
                continue;
            emplace_back(coeffs[i], roots[i], powers[i]);
            if (i + 1 < n && std::abs(roots[i] - std::conj(roots[i + 1])) < epsilon && std::abs(
                    coeffs[i] - std::conj(coeffs[i + 1])) < epsilon)
                ++i;
        }
    }

    explicit LaplaceExpression(std::vector<std::unique_ptr<Term<Type>>>&& terms, Type init_value = 0) :
        init_value(init_value),
        terms(std::move(terms)) {
    }

    ~LaplaceExpression() = default;

    void emplace_back(const Comp& c, const Comp& r, int n) {
        if (std::abs(r) == 0)
            switch (n) {
                case 0:
                    init_value += c.real();
                    return;
                case 1:
                    terms.push_back(std::make_unique<TimeTerm<Type>>(c.real()));
                    return;
                default:
                    terms.push_back(std::make_unique<PolyTerm<Type>>(c.real(), n));
                    return;
            }

        if (r.imag() == 0)
            switch (n) {
                case 0:
                    terms.push_back(std::make_unique<ExpTerm<Type>>(c.real(), r.real()));
                    return;
                case 1:
                    terms.push_back(std::make_unique<ExpTimeTerm<Type>>(c.real(), r.real()));
                    return;
                default:
                    terms.push_back(std::make_unique<ExpPolyTerm<Type>>(c.real(), r.real(), n));
                    return;
            }

        if (r.real() == 0)
            switch (n) {
                case 0:
                    terms.push_back(std::make_unique<CosTerm<Type>>(c, r));
                    return;
                case 1:
                    terms.push_back(std::make_unique<CosTimeTerm<Type>>(c, r));
                    return;
                default:
                    terms.push_back(std::make_unique<CosPolyTerm<Type>>(c, r, n));
                    return;
            }

        switch (n) {
            case 0:
                terms.push_back(std::make_unique<ExpCosTerm<Type>>(c, r));
                return;
            case 1:
                terms.push_back(std::make_unique<ExpCosTimeTerm<Type>>(c, r));
                return;
            default:
                terms.push_back(std::make_unique<ExpCosPolyTerm<Type>>(c, r, n));
        }
    }

    void push_back(std::vector<std::unique_ptr<Term<Type>>>&& new_terms) {
        terms.insert(
            terms.end(),
            std::make_move_iterator(new_terms.begin()),
            std::make_move_iterator(new_terms.end())
            );
    }

    Type operator()(double t) const {
        Type result = init_value;
        for (const auto& term : terms)
            result += term->value(t);
        return result;
    }

    LaplaceExpression derivative() const {
        Type derivative_init_value = 0;
        std::vector<std::unique_ptr<Term<Type>>> derivative_terms;
        derivative_terms.reserve(terms.size() * 2);
        for (const auto& term : terms) {
            derivative_init_value += term->derivativeConstant();
            auto derivative       = term->derivative();
            derivative_terms.insert(
                derivative_terms.end(),
                std::make_move_iterator(derivative.begin()),
                std::make_move_iterator(derivative.end())
                );
        }
        return LaplaceExpression(derivative_terms, derivative_init_value);
    }

    [[nodiscard]] std::string string() const {
        if (init_value != 0) {
            std::string result = (std::stringstream() << init_value).str();
            for (const auto& term : terms)
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
        if (this != &other) {
            init_value = other.init_value;
            terms.clear();
            terms.reserve(other.terms.size());
            for (const auto& term : other.terms)
                terms.push_back(term->clone());
        }
        return *this;
    }

    LaplaceExpression& operator=(LaplaceExpression&& other) noexcept {
        if (this != &other) {
            init_value = other.init_value;
            terms      = std::move(other.terms);
        }
        return *this;
    }
};
}
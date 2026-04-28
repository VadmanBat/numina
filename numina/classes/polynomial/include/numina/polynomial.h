#pragma once
// Created by Vadim on 18.06.2025.
#include <vector>
#include <complex>

class Polynomial {
public:
    using Type        = double;
    using Complex     = std::complex<Type>;
    using LongType    = long double;
    using LongComplex = std::complex<LongType>;
    using MultiRoots  = std::pair<
        std::vector<std::pair<Type, std::size_t>>,
        std::vector<std::pair<Complex, std::size_t>>
    >;

private:
    std::vector<Type> coeffs;
    std::size_t n;
    Type* c;

    void remove_leading_zeros();
    void zeroing();

public:
    constexpr Polynomial(const Type& constant = {});
    Polynomial(const std::vector<Type>& coefficients);
    Polynomial(std::vector<Type>&& coefficients) noexcept;
    Polynomial(const Polynomial& other);
    Polynomial(Polynomial&& other) noexcept;

    ~Polynomial() = default;

    [[nodiscard]] constexpr int degree() const noexcept { return static_cast<int>(n) - 1; }
    [[nodiscard]] constexpr const Type& operator[](const std::size_t index) const noexcept { return c[index]; }
    [[nodiscard]] constexpr Type& operator[](const std::size_t index) noexcept { return c[index]; }
    [[nodiscard]] constexpr const Type* data() const noexcept { return c; }
    [[nodiscard]] constexpr Type* data() noexcept { return c; }
    [[nodiscard]] std::vector<Type> vector() const { return coeffs; }
    [[nodiscard]] std::vector<Type> extractVector() && noexcept { return std::move(coeffs); }
    [[nodiscard]] constexpr bool isZero() const noexcept { return n == 1 && c[0] == 0; }

    [[nodiscard]] Polynomial derivative() const;
    [[nodiscard]] Polynomial derivative(std::size_t k) const;
    [[nodiscard]] Polynomial integral(const Type& constant = {}) const;
    void makeDerivative();
    void makeDerivative(std::size_t k);
    void makeIntegral(const Type& constant = {});

    [[nodiscard]] Polynomial deflate(const Type& root) const;
    [[nodiscard]] Polynomial deflateConjRoot(const Complex& root) const;
    void makeDeflate(const Type& root);
    void makeDeflateConjRoot(const Complex& root);

    [[nodiscard]] Polynomial monic() const;
    void makeMonic();

    [[nodiscard]] Polynomial compose(const Polynomial& other) const;

    [[nodiscard]] std::vector<Complex> computeRoots() const;
    [[nodiscard]] MultiRoots computeMultiRoots() const;

    [[nodiscard]] Type multiplicity(const Type& x, std::size_t m = 1) const;
    [[nodiscard]] Type multiplicity(const Complex& x, std::size_t m = 1) const;
    [[nodiscard]] std::size_t computeMultiplicity(const Type& x) const;
    [[nodiscard]] std::size_t computeMultiplicity(const Complex& x) const;

    Type operator()(const Type& x) const;
    Complex operator()(const Complex& x) const;
    LongType operator()(const LongType& x) const;
    LongComplex operator()(const LongComplex& x) const;

    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;

    Polynomial& operator=(const std::vector<Type>& coefficients);
    Polynomial& operator=(std::vector<Type>&& coefficients) noexcept;
    Polynomial& operator=(const Polynomial& other);
    Polynomial& operator=(Polynomial&& other) noexcept;
    Polynomial& operator=(const Type& constant);

    Polynomial operator+() const;
    Polynomial operator-() const;
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;
    Polynomial operator/(const Polynomial& divisor) const;
    Polynomial operator%(const Polynomial& divisor) const;

    [[nodiscard]] std::pair<Polynomial, Polynomial> divmod(const Polynomial& divisor) const;

    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);
    Polynomial& operator/=(const Polynomial& divisor);
    Polynomial& operator%=(const Polynomial& divisor);

    Polynomial& operator+=(const Type& scalar);
    Polynomial& operator-=(const Type& scalar);
    Polynomial& operator*=(const Type& scalar);
    Polynomial& operator/=(const Type& scalar);

    friend Polynomial operator+(const Polynomial& poly, const Type& scalar);
    friend Polynomial operator-(const Polynomial& poly, const Type& scalar);
    friend Polynomial operator*(const Polynomial& poly, const Type& scalar);
    friend Polynomial operator/(const Polynomial& poly, const Type& scalar);

    friend Polynomial operator-(const Type& scalar, const Polynomial& poly);
    friend Polynomial operator+(const Type& scalar, const Polynomial& poly);
    friend Polynomial operator*(const Type& scalar, const Polynomial& poly);

    friend Polynomial gcd(Polynomial a, Polynomial b);

    friend std::ostream& operator<<(std::ostream& os, const Polynomial& poly);
    friend std::istream& operator>>(std::istream& is, Polynomial& poly);
};

constexpr Polynomial::Polynomial(const Type& constant) :
    coeffs(1, constant),
    n(1),
    c(coeffs.data()) {
}
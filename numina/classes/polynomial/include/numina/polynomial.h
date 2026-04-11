//
// Created by Vadim on 18.06.2025.
//

#ifndef NUMINA_POLYNOMIAL_H
#define NUMINA_POLYNOMIAL_H

#include <vector>
#include <complex>

//template <typename Type>
class Polynomial {
private:
    using Type = double;
    using Complex = std::complex <Type>;
    using Roots = std::pair <std::vector <std::pair <Type, int>>, std::vector <std::pair <Complex, int>>>;

    std::vector <Type> coeffs;
    std::size_t n;
    Type* c;

    void remove_leading_zeros();
    void zeroing();

public:
    Polynomial(const Type& constant = Type(0));
    Polynomial(const std::vector <Type>& coefficients);
    Polynomial(std::vector <Type>&& coefficients) noexcept;
    Polynomial(const Polynomial& other);
    Polynomial(Polynomial&& other) noexcept;

    ~Polynomial() = default;

    [[nodiscard]] inline int degree() const { return static_cast<int>(n) - 1; }
    [[nodiscard]] inline const Type& operator[](std::size_t index) const { return c[index]; }
    [[nodiscard]] inline Type& operator[](std::size_t index) { return c[index]; }
    [[nodiscard]] inline const Type* data() const { return c; }
    [[nodiscard]] inline Type* data() { return c; }
    [[nodiscard]] inline std::vector <Type> vector() const { return coeffs; }
    [[nodiscard]] inline bool isZero() const { return n == 1 && c[0] == 0; }

    [[nodiscard]] Polynomial derivative() const;
    [[nodiscard]] Polynomial derivative(int k) const;
    [[nodiscard]] Polynomial integral(const Type& constant = Type(0)) const;
    void makeDerivative();
    void makeDerivative(int k);
    void makeIntegral(const Type& constant = Type(0));

    [[nodiscard]] Polynomial deflate(const Type& root) const;
    [[nodiscard]] Polynomial deflateConjRoot(const Complex& root) const;
    void makeDeflate(const Type& root);
    void makeDeflateConjRoot(const Complex& root);

    [[nodiscard]] Polynomial monic() const;
    void makeMonic();

    [[nodiscard]] Polynomial compose(const Polynomial& other) const;

    [[nodiscard]] std::vector <Complex> computeRoots() const;
    [[nodiscard]] Roots computeRootsWithMultiplicity() const;
    [[nodiscard]] bool isRoot(const Type& root, const Type& tolerance = 1e-6) const;
    [[nodiscard]] bool isRoot(const Complex& root, const Type& tolerance = 1e-6) const;
    [[nodiscard]] Type refineRoot(Type root, const Type& tolerance = 1e-6) const;
    [[nodiscard]] Complex refineRoot(Complex root, const Type& tolerance = 1e-6) const;

    [[nodiscard]] Type multiplicity(const Type& x) const;
    [[nodiscard]] Type multiplicity(const Complex& x) const;
    [[nodiscard]] int computeMultiplicity(const Type& x) const;
    [[nodiscard]] int computeMultiplicity(const Complex& x) const;

    Type operator()(const Type& x) const;
    Complex operator()(const Complex& x) const;
    bool operator==(const Polynomial& other) const;
    bool operator!=(const Polynomial& other) const;

    Polynomial& operator=(const std::vector <Type>& coefficients);
    Polynomial& operator=(std::vector <Type>&& coefficients) noexcept;
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

    [[nodiscard]] std::pair <Polynomial, Polynomial> divmod(const Polynomial& divisor) const;

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

#endif //NUMINA_POLYNOMIAL_H

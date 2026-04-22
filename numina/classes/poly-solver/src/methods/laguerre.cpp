#include "numina/poly-solver.h"
// Created by Vadim on 21.08.2025.
namespace numina {
std::vector<PolySolver::Complex> PolySolver::solve() {
    std::vector<Complex> roots;
    roots.reserve(degree);

    Type error_previous, error_current;
    Complex x, fx, d1fx, d2fx, b, discriminant, a1, a2, gx, root;

    while (degree > 0) {
        x = Complex(-1, -1);
        for (int i = 0; i < LAGUERRE_MAX_ITER; ++i) {
            fx = f(x);
            if (std::abs(fx) < E1 * std::abs(x))
                break;

            d1fx         = d1(x);
            d2fx         = d2(x);
            b            = static_cast<Type>(degree - 1) * d1fx;
            discriminant = std::sqrt(b * b - static_cast<Type>(degree * (degree - 1)) * fx * d2fx);
            a1           = d1fx - discriminant;
            a2           = d1fx + discriminant;
            x            -= static_cast<Type>(degree) * fx / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        }

        int m = 1;
        if (std::abs(fx) < E2 * std::abs(x)) {
            m = df[0].computeMultiplicity(x);
            if (m > 2)
                for (std::size_t i = df.size() - 1; i < m; ++i)
                    df.emplace_back(df[i].derivative());
        }

        const auto& g  = df[m - 1];
        const auto& dg = df[m];

        root           = x;
        gx             = g(x);
        error_previous = std::abs(gx);
        for (int i = 0; i < NEWTON_MAX_ITER; ++i) {
            x             -= gx / dg(x);
            gx            = g(x);
            error_current = std::abs(gx);
            if (error_current >= error_previous) {
                if (error_current == error_previous)
                    root = x;
                break;
            }
            error_previous = error_current;
            root           = x;
        }

        if (std::abs(root.imag()) > E1 * std::abs(root)) {
            for (int i = 0; i < m; ++i) {
                roots.emplace_back(root);
                roots.emplace_back(std::conj(root));
            }
            deflate_conj(root, m);
        }
        else {
            const auto r = root.real();
            roots.resize(roots.size() + m, r);
            deflate(r, m);
        }
    }

    clear();
    return roots;
}

void PolySolver::solve_with_multiplicities() {
    Type error_previous, error_current;
    Complex x, fx, d1fx, d2fx, b, discriminant, a1, a2, gx, root;

    while (degree > 0) {
        x = Complex(-1, -1);
        for (int i = 0; i < LAGUERRE_MAX_ITER; ++i) {
            fx = f(x);
            if (std::abs(fx) < E1 * std::abs(x))
                break;

            d1fx         = d1(x);
            d2fx         = d2(x);
            b            = static_cast<Type>(degree - 1) * d1fx;
            discriminant = std::sqrt(b * b - static_cast<Type>(degree * (degree - 1)) * fx * d2fx);
            a1           = d1fx - discriminant;
            a2           = d1fx + discriminant;
            x            -= static_cast<Type>(degree) * fx / (std::abs(a1) > std::abs(a2) ? a1 : a2);
        }

        int m = 1;
        if (std::abs(fx) < E2 * std::abs(x)) {
            m = df[0].computeMultiplicity(x);
            if (m > 2)
                for (std::size_t i = df.size() - 1; i < m; ++i)
                    df.emplace_back(df[i].derivative());
        }

        const auto& g  = df[m - 1];
        const auto& dg = df[m];

        root           = x;
        gx             = g(x);
        error_previous = std::abs(gx);
        for (int i = 0; i < NEWTON_MAX_ITER; ++i) {
            x             -= gx / dg(x);
            gx            = g(x);
            error_current = std::abs(gx);
            if (error_current >= error_previous) {
                if (error_current == error_previous)
                    root = x;
                break;
            }
            error_previous = error_current;
            root           = x;
        }

        if (std::abs(root.imag()) > E1 * std::abs(root)) {
            result.second.emplace_back(root, m);
            deflate_conj(root, m);
        }
        else {
            const auto r = root.real();
            result.first.emplace_back(r, m);
            deflate(r, m);
        }
    }
}

void PolySolver::solve_with_implicit_deflation() {
    std::vector<std::pair<Complex, int>> found; // (корень, кратность)
    std::size_t m_eff = degree;

    Complex x, fx, d1fx, d2fx, G0, H0, S, T, G, H, inside, sd, denom1, denom2, denom;

    while (m_eff > 0) {
        // ==============================================
        // 1. ЛАГЕРР 4-го ПОРЯДКА С НЕЯВНОЙ ДЕФЛЯЦИЕЙ
        // ==============================================
        x = Complex(-1.0, -1.0); // стартовая точка (можно улучшить позже)

        for (int i = 0; i < LAGUERRE_MAX_ITER; ++i) {
            fx = f(x);
            if (std::abs(fx) < E1 * std::abs(x))
                break;

            d1fx = d1(x);
            d2fx = d2(x);

            // === НЕЯВНАЯ ДЕФЛЯЦИЯ (implicit deflation) ===
            S = T = Complex(0.0, 0.0);
            for (const auto& [r, mu] : found) {
                Complex den = x - r;
                if (std::abs(den) < E1)
                    den = Complex(E1, 0.0);
                Complex inv = Type(1.0) / den;
                S += static_cast<Type>(mu) * inv;
                T += static_cast<Type>(mu) * inv * inv;
            }

            G0 = d1fx / fx;
            H0 = G0 * G0 - d2fx / fx;
            G  = G0 - S;
            H  = H0 - T;

            // === ФОРМУЛА ЛАГЕРРА 4-го ПОРЯДКА ===
            inside = static_cast<Type>(m_eff) * H - G * G;
            inside = static_cast<Type>(m_eff - 1) * inside;
            sd     = std::sqrt(inside);

            denom1 = G + sd;
            denom2 = G - sd;
            denom  = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;

            if (std::abs(denom) < E2)
                denom = Complex(E2, 0.0);

            x -= static_cast<Type>(m_eff) / denom;
        }

        // ==============================================
        // 2. Поиск кратности
        // ==============================================
        int m = 1;
        if (std::abs(f(x)) < E2 * std::abs(x)) {
            m = df[0].computeMultiplicity(x);
            if (m > 2)
                for (std::size_t i = df.size() - 1; i < static_cast<std::size_t>(m); ++i)
                    df.emplace_back(df[i].derivative());
        }

        // ==============================================
        // 3. Уточнение (МОДИФИЦИРОВАННЫЙ НЬЮТОН С implicit)
        // ==============================================
        Complex root = x;
        const auto& g  = df[m - 1];
        const auto& dg = df[m];

        Complex gx = g(x);

        for (int i = 0; i < NEWTON_MAX_ITER; ++i) {
            Complex dg_val = dg(x);

            // НЕЯВНАЯ ДЕФЛЯЦИЯ в уточнении
            Complex S(0.0);
            for (const auto& [r, mu] : found) {
                if (mu < m) continue;
                Complex den = x - r;
                S += static_cast<Type>(mu - m + 1) / den;
            }

            Complex denom = dg_val - gx * S;          // ← implicit
            if (std::abs(denom) < E2) break;

            Complex step = gx / denom;                // ← используем implicit!
            x -= step;
            gx = g(x);

            if (std::abs(step) < E1 * (Type(1.0) + std::abs(x))) break;
        }

        root = x;

        // ==============================================
        // 4. Сохранение корня (implicit + сопряжённый)
        // ==============================================
        found.emplace_back(root, m);
        m_eff -= static_cast<std::size_t>(m);

        if (std::abs(root.imag()) > E1 * std::abs(root)) {
            Complex conj_root = std::conj(root);
            found.emplace_back(conj_root, m);
            m_eff -= static_cast<std::size_t>(m);

            result.second.emplace_back(root, m);
            result.second.emplace_back(conj_root, m);
        }
        else {
            result.first.emplace_back(root.real(), m);
        }
    }
}
/*
PolySolver::Roots PolySolver::solve_with_implicit_deflation() {
    Roots roots;

    if (degree == 0) {
        clear_state();
        return roots;
    }
    if (degree == 1) {
        Type root = -coeffs[1] / coeffs[0];
        roots.first.emplace_back(root, 1);
        clear_state();
        return roots;
    }

    std::vector<std::pair<Complex, int>> found; // (корень, кратность)
    std::size_t m_eff = degree;

    size_t zero_mult = 0;
    while (zero_mult < degree && coeffs[degree - zero_mult] == 0)
        ++zero_mult;

    if (zero_mult > 0) {
        found.emplace_back(Complex(0.0, 0.0), static_cast<int>(zero_mult));
        m_eff -= zero_mult;

        // если весь многочлен был x^k — сразу возвращаем
        if (m_eff == 0) {
            roots.first.emplace_back(0.0, static_cast<int>(zero_mult));
            clear_state();
            return roots;
        }
    }

    Complex x, fx, d1fx, d2fx, G0, H0, S, T, G, H, inside, sd, denom1, denom2, denom;

    while (m_eff > 0) {
        //std::cout << "loop-1\n";
        // ==============================================
        // 1. Лаггер с НЕЯВНОЙ дефляцией (ИСПРАВЛЕННАЯ формула)
        // ==============================================
        x = Complex(-1.0, -1.0); // можно потом улучшить, как в оригинальном solve()
        for (int i = 0; i < 30; ++i) {
            fx = f(x);
            if (std::abs(fx) < e1 * std::abs(x))
                break;

            d1fx = d1(x);
            d2fx = d2(x);

            // === НЕЯВНАЯ ДЕФЛЯЦИЯ (всё так же) ===
            S = T = Complex(0.0);
            for (const auto& [r, mu] : found) {
                Complex den = x - r;
                if (std::abs(den) < e1)
                    den = Complex(e1, 0.0);
                Complex inv = Type(1.0) / den;
                S           += static_cast<Type>(mu) * inv;
                T           += static_cast<Type>(mu) * inv * inv;
            }

            G0 = d1fx / fx;
            H0 = G0 * G0 - d2fx / fx;
            G  = G0 - S;
            H  = H0 - T;

            // === ПРАВИЛЬНАЯ ФОРМУЛА ЛАГЕРРА ===
            inside = static_cast<Type>(m_eff) * H - G * G;
            inside = static_cast<Type>(m_eff - 1) * inside; // ← вот здесь было главное отличие
            sd     = std::sqrt(inside);

            denom1 = G + sd;
            denom2 = G - sd;
            denom  = (std::abs(denom1) > std::abs(denom2)) ? denom1 : denom2;

            if (std::abs(denom) < e2)
                denom = Complex(e2, 0.0);

            x -= static_cast<Type>(m_eff) / denom;
        }

        // ==============================================
        // 2. Поиск кратности (без изменений)
        // ==============================================
        int m = 1;
        if (std::abs(f(x)) < e2 * std::abs(x)) {
            m = compute_multiplicity(x);
            if (m > 2)
                for (std::size_t i = dfx.size() - 1; i < static_cast<std::size_t>(m); ++i)
                    dfx.emplace_back(dfx[i].derivative());
        }

        // ==============================================
        // 3. Уточнение (без изменений)
        // ==============================================
        Complex root = x;
        // g = p^{(m-1)}, dg = p^{(m)}
        const auto& g  = dfx[m - 1];
        const auto& dg = dfx[m];

        Complex gx = g(x);

        for (int i = 0; i < 30; ++i) {   // больше итераций
            Complex dg_val = dg(x);

            // === МОДИФИЦИРОВАННЫЙ НЬЮТОН С НЕЯВНОЙ ДЕФЛЯЦИЕЙ ===
            Complex S(0.0);
            for (const auto& [r, mu] : found) {
                if (mu < m)
                    continue;

                Complex den = x - r;
                //if (std::abs(den) < e1) den = Complex(e1, 0.0);
                S += static_cast<Type>(mu - m + 1) / den;
            }

            Complex denom = dg_val - gx * S;   // ← ключевой момент
            //if (std::abs(denom) < e2) break;

            Complex step = gx / denom;
            step = gx / dg_val;
            x -= step;
            gx = g(x);

            //if (std::abs(step) < e1 * (Type(1.0) + std::abs(x))) break;
        }

        root = x;

        // ==============================================
        // 4. Запоминаем корень (НЕЯВНАЯ дефляция + сопряжённый)
        // ==============================================
        found.emplace_back(root, m);
        m_eff -= static_cast<std::size_t>(m);
        //std::cout << "m_eff = " << m_eff << std::endl;

        if (std::abs(root.imag()) > e1 * std::abs(root)) {
            // Комплексный корень — сразу добавляем сопряжённый
            Complex conj_root = std::conj(root);
            found.emplace_back(conj_root, m);
            m_eff -= static_cast<std::size_t>(m);

            roots.second.emplace_back(root, m);
            roots.second.emplace_back(conj_root, m);
        }
        else
            roots.first.emplace_back(root.real(), m);
    }

    clear_state();
    return roots;
}*/
}
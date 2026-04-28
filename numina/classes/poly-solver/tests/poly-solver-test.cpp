#include "numina/poly-solver.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

using Type    = double;
using Complex = std::complex<Type>;

// Horner evaluation (high-to-low coefficients)
Type evaluate_poly(const std::vector<Type>& coeffs, Complex z) {
    if (coeffs.empty()) return 0.0;
    Complex result = coeffs[0];
    for (size_t i = 1; i < coeffs.size(); ++i)
        result = result * z + coeffs[i];
    return std::abs(result);
}

// Compute max absolute coefficient (for relative error)
Type max_abs_coeff(const std::vector<Type>& coeffs) {
    Type max_val = 0.0;
    for (auto c : coeffs)
        max_val = std::max(max_val, std::abs(c));
    return max_val > 0 ? max_val : 1.0;
}

// Print one test result
void print_test_result(int degree,
                       double max_err,
                       double avg_err,
                       double rel_err,
                       int root_count,
                       bool has_nan,
                       bool has_inf,
                       bool print_roots,
                       const std::vector<std::pair<Type, std::size_t>>& real_roots,
                       const std::vector<std::pair<Complex, std::size_t>>& complex_roots) {
    std::cout << std::left << std::setw(6) << ("deg-" + std::to_string(degree))
              << "  Max error: " << std::scientific << std::setprecision(2) << max_err
              << "   Avg error: " << std::scientific << std::setprecision(2) << avg_err
              << "   Rel error: " << std::scientific << std::setprecision(2) << rel_err
              << "   root-count = " << root_count;

    if (has_nan || has_inf) {
        std::cout << "   [";
        if (has_nan) std::cout << "NaN";
        if (has_nan && has_inf) std::cout << "/";
        if (has_inf) std::cout << "Inf";
        std::cout << " detected!]";
    }
    std::cout << std::endl;

    if (print_roots) {
        for (const auto& [r, mult] : real_roots) {
            std::cout << "  " << r << " (" << mult << ")\n";
        }
        for (const auto& [r, mult] : complex_roots) {
            std::cout << "  " << r << " (" << mult << ")\n";
        }
        std::cout << '\n';
    }
}

// Returns pair<has_nan, has_inf> for the whole group
std::pair<bool, bool> run_group(const std::string& filename, numina::PolySolver& solver, bool print_roots = false) {
    std::string path = "tests/data/" + filename;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open " << path << std::endl;
        return {false, false};
    }

    std::string group_name = filename.substr(0, filename.find_last_of('.'));
    std::cout << "\n=== " << group_name << " ===\n" << std::string(65, '-') << std::endl;

    int N;
    file >> N;

    bool group_has_nan = false;
    bool group_has_inf = false;

    for (int test_idx = 0; test_idx < N; ++test_idx) {
        int n;
        file >> n;
        std::vector<Type> coeffs(n + 1);
        for (int i = 0; i <= n; ++i) {
            file >> coeffs[i];
        }

        auto roots = solver.multisolve(coeffs, numina::PolySolver::Method::Implicit);

        double max_err = 0.0;
        double sum_err = 0.0;
        int count = 0;
        int root_count = 0;
        bool has_nan = false;
        bool has_inf = false;

        for (const auto& [r, mult] : roots.first) {
            double err = evaluate_poly(coeffs, r);
            if (std::isnan(err)) has_nan = true;
            if (std::isinf(err)) has_inf = true;
            max_err = std::max(max_err, err);
            sum_err += err;
            count++;
            root_count += mult;
        }
        for (const auto& [r, mult] : roots.second) {
            double err = evaluate_poly(coeffs, r);
            if (std::isnan(err)) has_nan = true;
            if (std::isinf(err)) has_inf = true;
            max_err = std::max(max_err, err);
            sum_err += err;
            count++;
            root_count += mult;
        }

        if (has_nan) group_has_nan = true;
        if (has_inf) group_has_inf = true;

        double avg_err = (count > 0) ? sum_err / count : 0.0;
        double max_coeff = max_abs_coeff(coeffs);
        double rel_err = max_err / max_coeff;

        print_test_result(n, max_err, avg_err, rel_err, root_count,
                          has_nan, has_inf, print_roots, roots.first, roots.second);
    }

    return {group_has_nan, group_has_inf};
}

int main(int argc, char* argv[]) {
    numina::PolySolver solver;

    bool print_roots = false;
    if (argc > 1 && std::string(argv[1]) == "--print-roots") {
        print_roots = true;
    }

    std::cout << "=== PolySolver Implicit Deflation - File Based Tests ===\n";
    std::cout << "Testing only the specified 10 groups from tests/data/\n\n";

    const std::vector<std::string> test_files = {
        "random.txt",
        "complex-roots.txt",
        "wilkinson.txt",
        "multiple-roots.txt",
        "close-roots.txt",
        "close-cluster.txt",
        "small-coeffs.txt",
        "large-scale.txt",
        "large-magnitude.txt",
        "extreme-coeffs.txt",
    };

    std::vector<std::string> nan_groups;
    std::vector<std::string> inf_groups;

    for (const auto& f : test_files) {
        auto [has_nan, has_inf] = run_group(f, solver, print_roots);
        if (has_nan) nan_groups.push_back(f);
        if (has_inf) inf_groups.push_back(f);
    }

    std::cout << "\n" << std::string(75, '=') << std::endl;
    std::cout << "All tests completed.\n\n";

    // === NaN / Inf Summary ===
    std::cout << "=== NaN / Inf Summary ===\n";
    if (!nan_groups.empty()) {
        std::cout << "NaN detected in: ";
        for (size_t i = 0; i < nan_groups.size(); ++i) {
            std::cout << nan_groups[i];
            if (i + 1 < nan_groups.size()) std::cout << ", ";
        }
        std::cout << "\n";
    } else {
        std::cout << "No NaN detected in any test group.\n";
    }

    if (!inf_groups.empty()) {
        std::cout << "Inf detected in: ";
        for (size_t i = 0; i < inf_groups.size(); ++i) {
            std::cout << inf_groups[i];
            if (i + 1 < inf_groups.size()) std::cout << ", ";
        }
        std::cout << "\n";
    } else {
        std::cout << "No Inf detected in any test group.\n";
    }

    return 0;
}
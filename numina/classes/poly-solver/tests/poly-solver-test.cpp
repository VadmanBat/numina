#include "numina/poly-solver.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

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

void run_group(const std::string& filename, numina::PolySolver& solver) {
    std::string path = "tests/data/" + filename;
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cout << "ERROR: Cannot open " << path << std::endl;
        return;
    }

    std::string group_name = filename.substr(0, filename.find_last_of('.'));
    std::cout << "\n=== " << group_name << " ===\n" << std::string(65, '-') << std::endl;

    int N;
    file >> N;

    for (int test_idx = 0; test_idx < N; ++test_idx) {
        int n;
        file >> n;
        std::vector<Type> coeffs(n + 1);
        for (int i = 0; i <= n; ++i) {
            file >> coeffs[i];
        }

        auto roots = solver.solveWithMultiplicities(coeffs);

        double max_err = 0.0;
        double sum_err = 0.0;
        int count = 0;
        int root_count = 0;

        for (const auto& [r, mult] : roots.first) {
            double err = evaluate_poly(coeffs, r);
            max_err = std::max(max_err, err);
            sum_err += err;
            count++;
            root_count += mult;
        }
        for (const auto& [r, mult] : roots.second) {
            double err = evaluate_poly(coeffs, r);
            max_err = std::max(max_err, err);
            sum_err += err;
            count++;
            root_count += mult;
        }

        double avg_err = (count > 0) ? sum_err / count : 0.0;

        std::cout << std::left << std::setw(6) << ("deg-" + std::to_string(n))
                  << "  Max error: " << std::scientific << std::setprecision(2) << max_err
                  << "   Avg error: " << std::scientific << std::setprecision(2) << avg_err
                  << "   root-count = " << root_count
                  << std::endl;

        /*for (const auto& [r, mult] : roots.first) {
            std::cout << r << " (" << mult << ")\n";
        }
        for (const auto& [r, mult] : roots.second) {
            std::cout << r << " (" << mult << ")\n";
        }
        std::cout << '\n';*/
    }
}

int main() {
    numina::PolySolver solver;

    std::cout << "=== PolySolver Implicit Deflation - File Based Tests ===\n";
    std::cout << "Testing only the specified 10 groups from tests/data/\n\n";

    std::vector<std::string> test_files = {
        "random.txt",
        "complex-roots.txt",
        "wilkinson.txt",
        "multiple-roots.txt",
        "close-roots.txt",
        "close-cluster.txt",
        "small-coeffs.txt",
        "large-scale.txt",
        //"large-magnitude.txt",
        //"extreme-coeffs.txt",
    };

    for (const auto& f : test_files) {
        run_group(f, solver);
    }

    std::cout << "\n" << std::string(75, '=') << std::endl;
    std::cout << "All tests completed.\n";
    std::cout << "Good performance: Max error < 1e-8, Avg error < 1e-10\n";

    return 0;
}
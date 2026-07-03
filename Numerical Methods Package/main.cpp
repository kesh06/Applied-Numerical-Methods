#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "linear.hpp"
#include "interpolation.hpp"
#include "non_linear.hpp"
#include "quadrature.hpp"

// =====================================================================
// Global Mathematical Test Helpers
// =====================================================================
double f_nonlinear(double x) { return x * x - 2.0; }          // Root at sqrt(2) ≈ 1.414214
double df_nonlinear(double x) { return 2.0 * x; }
double f_flat(double x) { return (x - 2.0) * (x - 2.0) + 1.0; } // No real roots (flat turning point near 2)
double f_quadrature(double x) { return std::sin(x); }          // Integral from 0 to pi is exactly 2.0

void printVector(const std::vector<double>& vec) {
    std::cout << "[ ";
    for (double val : vec) std::cout << val << " ";
    std::cout << "]\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "==================================================================\n";
    std::cout << "                NUMERICAL SUITE VALIDATION BENCHMARK              \n";
    std::cout << "==================================================================\n\n";

    // =====================================================================
    // SECTION 1: LINEAR SYSTEM MODULE
    // =====================================================================
    std::cout << "--------------------------------------------------\n";
    std::cout << "1. TESTING: LINEAR MODULE (Direct, Iterative, Norms)\n";
    std::cout << "--------------------------------------------------\n";
    
    std::vector<std::vector<double>> A_valid = {
        {4.0, -1.0, 1.0},
        {-1.0, 4.0, -2.0},
        {1.0, -2.0, 4.0}
    };
    std::vector<double> b_valid = {12.0, -1.0, 5.0}; // Expected solution: [3, 1, 1]

    std::cout << "Matrix Inf-Norm: " << NumericalMethods::Linear::norm_inf(A_valid) << "\n";
    std::cout << "Matrix L1-Norm:  " << NumericalMethods::Linear::norm_l1(A_valid) << "\n";
    std::cout << "Matrix L2-Norm:  " << NumericalMethods::Linear::norm_l2(A_valid) << "\n\n";

    try {
        auto x_gauss = NumericalMethods::Linear::gauss_scaled(A_valid, b_valid);
        std::cout << "Gauss Scaled Solution:   "; printVector(x_gauss);
        
        auto x_sor = NumericalMethods::Linear::sor(A_valid, b_valid, 1.15, 1e-7);
        std::cout << "SOR Iterative Solution:  "; printVector(x_sor);
    } catch (const std::exception& e) {
        std::cout << "Unexpected Linear Solver Failure: " << e.what() << "\n";
    }

    std::cout << "\nTesting Singular Matrix Boundary Case:\n";
    std::vector<std::vector<double>> A_singular = {{1, 2}, {2, 4}};
    std::vector<double> b_singular = {5, 10};
    try {
        auto x_fail = NumericalMethods::Linear::gauss(A_singular, b_singular);
        std::cout << "Failure: System allowed computation on a singular matrix.\n";
    } catch (const std::exception& e) {
        std::cout << "Caught Expected Exception: " << e.what() << "\n";
    }
    std::cout << "\n";

    // =====================================================================
    // SECTION 2: INTERPOLATION MODULE
    // =====================================================================
    std::cout << "--------------------------------------------------\n";
    std::cout << "2. TESTING: INTERPOLATION MODULE\n";
    std::cout << "--------------------------------------------------\n";
    
    std::vector<double> x_nodes = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> y_nodes = {2.0, 3.0, 5.0, 4.0};
    double target_x = 2.5;

    try {
        double y_newton = NumericalMethods::Interpolation::newton(x_nodes, y_nodes, target_x);
        double y_lagrange = NumericalMethods::Interpolation::lagrange(x_nodes, y_nodes, target_x);
        double y_spline = NumericalMethods::Interpolation::cubic_spline(x_nodes, y_nodes, target_x);

        std::cout << "Evaluating interpolants at target x = " << target_x << "\n";
        std::cout << "Newton Divided Difference: " << y_newton << "\n";
        std::cout << "Lagrange Polynomial:       " << y_lagrange << "\n";
        std::cout << "Natural Cubic Spline:      " << y_spline << "\n";
    } catch (const std::exception& e) {
        std::cout << "Unexpected Interpolation Failure: " << e.what() << "\n";
    }
    std::cout << "\n";

    // =====================================================================
    // SECTION 3: NON-LINEAR SYSTEM MODULE (ROOT FINDING)
    // =====================================================================
    std::cout << "--------------------------------------------------\n";
    std::cout << "3. TESTING: NON-LINEAR ROOT MODULE\n";
    std::cout << "--------------------------------------------------\n";
    
    try {
        double root_bis = NumericalMethods::NonLinear::bisection(f_nonlinear, 1.0, 2.0, 1e-6);
        double root_sec = NumericalMethods::NonLinear::secant(f_nonlinear, 1.0, 2.0, 1e-6);
        double root_nr  = NumericalMethods::NonLinear::newton_raphson(f_nonlinear, df_nonlinear, 2.0, 1e-6);

        std::cout << "Target function f(x) = x^2 - 2 (Exact = 1.414214):\n";
        std::cout << "Bisection Method Root:      " << root_bis << "\n";
        std::cout << "Secant Method Root:         " << root_sec << "\n";
        std::cout << "Newton-Raphson Method Root: " << root_nr << "\n";
    } catch (const std::exception& e) {
        std::cout << "Unexpected Root Finder Failure: " << e.what() << "\n";
    }

    std::cout << "\nTesting Flat Slope Exception Handling:\n";
    try {
        double root_flat = NumericalMethods::NonLinear::newton_raphson(f_flat, df_nonlinear, 2.0, 1e-6);
        std::cout << "Failure: System allowed zero-derivative step computation.\n";
    } catch (const std::exception& e) {
        std::cout << "Caught Expected Exception: " << e.what() << "\n";
    }
    std::cout << "\n";

    // =====================================================================
    // SECTION 4: QUADRATURE MODULE (INTEGRATION)
    // =====================================================================
    std::cout << "--------------------------------------------------\n";
    std::cout << "4. TESTING: QUADRATURE MODULE\n";
    std::cout << "--------------------------------------------------\n";
    
    double pi = std::acos(-1.0);
    int intervals = 100;

    try {
        double int_trap = NumericalMethods::Quadrature::trapezoidal(f_quadrature, 0.0, pi, intervals);
        double int_simp = NumericalMethods::Quadrature::simpson(f_quadrature, 0.0, pi, intervals);

        std::cout << "Integrating sin(x) from 0 to pi (Exact = 2.000000):\n";
        std::cout << "Composite Trapezoidal (n=" << intervals << "): " << int_trap << "\n";
        std::cout << "Composite Simpson 1/3 (n=" << intervals << "): " << int_simp << "\n";
    } catch (const std::exception& e) {
        std::cout << "Unexpected Quadrature Failure: " << e.what() << "\n";
    }

    std::cout << "\nTesting Simpson Interval Constraints:\n";
    try {
        double int_odd = NumericalMethods::Quadrature::simpson(f_quadrature, 0.0, pi, 15);
        std::cout << "Failure: System evaluated Simpson's rule with odd interval count.\n";
    } catch (const std::exception& e) {
        std::cout << "Caught Expected Exception: " << e.what() << "\n";
    }

    std::cout << "==================================================================\n";
    std::cout << "                 ALL BENCHMARK SECTIONS COMPLETED                 \n";
    std::cout << "==================================================================\n";

    return 0;
}
#ifndef NUMERICAL_METHODS_NON_LINEAR_HPP
#define NUMERICAL_METHODS_NON_LINEAR_HPP

#include <functional>
#include <cmath>
#include <stdexcept>

namespace NumericalMethods {
    namespace NonLinear {

        // =====================================================================
        // 1. Bisection Method
        // =====================================================================
        inline double bisection(std::function<double(double)> f, double a, double b, double tol = 1e-6, int max_iter = 100) {
            double fa = f(a);
            double fb = f(b);
            
            if (fa * fb >= 0) {
                throw std::invalid_argument("Error: f(a) and f(b) must have opposite signs.");
            }

            double c = a;
            for (int iter = 0; iter < max_iter; iter++) {
                c = a + (b - a) / 2.0; // Avoids potential overflow compared to (a + b) / 2
                double fc = f(c);

                if (std::abs(fc) < 1e-15 || (b - a) / 2.0 < tol) {
                    return c;
                }

                if (fa * fc < 0) {
                    b = c;
                    fb = fc;
                } else {
                    a = c;
                    fa = fc;
                }
            }
            return c;
        }

        // =====================================================================
        // 2. Regula Falsi (False Position) Method
        // =====================================================================
        inline double regula_falsi(std::function<double(double)> f, double a, double b, double tol = 1e-6, int max_iter = 100) {
            double fa = f(a);
            double fb = f(b);

            if (fa * fb >= 0) {
                throw std::invalid_argument("Error: f(a) and f(b) must have opposite signs.");
            }

            double c = a;
            for (int iter = 0; iter < max_iter; iter++) {
                // Secant line intercept point on the x-axis
                c = (a * fb - b * fa) / (fb - fa);
                double fc = f(c);

                if (std::abs(fc) < tol) {
                    return c;
                }

                if (fa * fc < 0) {
                    b = c;
                    fb = fc;
                } else {
                    a = c;
                    fa = fc;
                }
            }
            return c;
        }

        // =====================================================================
        // 3. Secant Method
        // =====================================================================
        inline double secant(std::function<double(double)> f, double x0, double x1, double tol = 1e-6, int max_iter = 100) {
            double f0 = f(x0);
            double f1 = f(x1);

            for (int iter = 0; iter < max_iter; iter++) {
                if (std::abs(f1 - f0) < 1e-12) {
                    throw std::runtime_error("Error: Secant method division by zero (flat slope approximation).");
                }

                double x_next = x1 - f1 * (x1 - x0) / (f1 - f0);

                if (std::abs(x_next - x1) < tol) {
                    return x_next;
                }

                x0 = x1;
                f0 = f1;
                x1 = x_next;
                f1 = f(x_next);
            }
            return x1;
        }

        // =====================================================================
        // 4. Newton-Raphson Method
        // =====================================================================
        inline double newton_raphson(std::function<double(double)> f, std::function<double(double)> df, double x0, double tol = 1e-6, int max_iter = 100) {
            double x = x0;

            for (int iter = 0; iter < max_iter; iter++) {
                double y = f(x);
                double dy = df(x);

                if (std::abs(dy) < 1e-12) {
                    throw std::runtime_error("Error: Newton-Raphson derivative approached zero (local extremum).");
                }

                double x_next = x - y / dy;

                if (std::abs(x_next - x) < tol) {
                    return x_next;
                }
                x = x_next;
            }
            return x;
        }
    }
}

#endif
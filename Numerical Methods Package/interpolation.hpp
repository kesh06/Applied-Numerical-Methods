#ifndef NUMERICAL_METHODS_INTERPOLATION_HPP
#define NUMERICAL_METHODS_INTERPOLATION_HPP

#include <vector>
#include <stdexcept>
#include <algorithm>
#include "linear.hpp" // Integrated to use the Gaussian solver for Splines

namespace NumericalMethods {
    namespace Interpolation {

        // =====================================================================
        // 1. Newton's Divided Difference Interpolation
        // =====================================================================
        inline double newton(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double x_eval) {
            int n = x_nodes.size();
            if (n == 0 || n != (int)y_nodes.size()) {
                throw std::invalid_argument("Error: Node vectors must be non-empty and of equal size.");
            }

            // Create divided difference table coefficients
            std::vector<double> coef = y_nodes;
            for (int j = 1; j < n; j++) {
                for (int i = n - 1; i >= j; i--) {
                    coef[i] = (coef[i] - coef[i - 1]) / (x_nodes[i] - x_nodes[i - j]);
                }
            }

            // Evaluate the polynomial using Horner's-like method
            double result = 0.0;
            for (int i = n - 1; i >= 0; i--) {
                result = coef[i] + (x_eval - x_nodes[i]) * result;
            }
            return result;
        }

        // =====================================================================
        // 2. Lagrange Interpolation
        // =====================================================================
        inline double lagrange(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double x_eval) {
            int n = x_nodes.size();
            if (n == 0 || n != (int)y_nodes.size()) {
                throw std::invalid_argument("Error: Node vectors must be non-empty and of equal size.");
            }

            double total_sum = 0.0;
            for (int i = 0; i < n; i++) {
                double term = y_nodes[i];
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        if (std::abs(x_nodes[i] - x_nodes[j]) < 1e-12) {
                            throw std::runtime_error("Error: Duplicate x-coordinates detected in nodes.");
                        }
                        term *= (x_eval - x_nodes[j]) / (x_nodes[i] - x_nodes[j]);
                    }
                }
                total_sum += term;
            }
            return total_sum;
        }

        // =====================================================================
        // 3. Natural Cubic Spline Interpolation
        // =====================================================================
        inline double cubic_spline(const std::vector<double>& x_nodes, const std::vector<double>& y_nodes, double x_eval) {
            int n = x_nodes.size() - 1; // Number of intervals
            if (n < 2 || x_nodes.size() != y_nodes.size()) {
                throw std::invalid_argument("Error: Spline requires at least 3 points.");
            }

            // Calculate step sizes h
            std::vector<double> h(n);
            for (int i = 0; i < n; i++) {
                h[i] = x_nodes[i + 1] - x_nodes[i];
                if (h[i] <= 0) {
                    throw std::invalid_argument("Error: x_nodes must be strictly increasing.");
                }
            }

            // Set up tridiagonal system system matrix A_sys and right side b_sys to find 'c' coefficients
            // The boundary conditions are natural: c[0] = 0 and c[n] = 0
            int num_eq = n + 1;
            std::vector<std::vector<double>> A_sys(num_eq, std::vector<double>(num_eq, 0.0));
            std::vector<double> b_sys(num_eq, 0.0);

            A_sys[0][0] = 1.0;
            A_sys[n][n] = 1.0;

            for (int i = 1; i < n; i++) {
                A_sys[i][i - 1] = h[i - 1];
                A_sys[i][i] = 2.0 * (h[i - 1] + h[i]);
                A_sys[i][i + 1] = h[i];
                b_sys[i] = 3.0 * (y_nodes[i + 1] - y_nodes[i]) / h[i] - 3.0 * (y_nodes[i] - y_nodes[i - 1]) / h[i - 1];
            }

            // Solve for c vector using the Scaled Partial Pivoting Gaussian routine from linear.hpp
            std::vector<double> c = NumericalMethods::Linear::gauss_scaled(A_sys, b_sys);

            // Compute remaining spline coefficients b and d for each interval
            std::vector<double> b(n), d(n);
            for (int i = 0; i < n; i++) {
                b[i] = (y_nodes[i + 1] - y_nodes[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
                d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
            }

            // Determine which interval x_eval falls inside (clipping boundaries safely)
            int idx = 0;
            if (x_eval <= x_nodes[0]) {
                idx = 0;
            } else if (x_eval >= x_nodes[n]) {
                idx = n - 1;
            } else {
                auto it = std::upper_bound(x_nodes.begin(), x_nodes.end(), x_eval);
                idx = std::distance(x_nodes.begin(), it) - 1;
            }

            // Evaluate the cubic polynomial for interval idx: S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3
            double dx = x_eval - x_nodes[idx];
            return y_nodes[idx] + b[idx] * dx + c[idx] * dx * dx + d[idx] * dx * dx * dx;
        }
    }
}

#endif
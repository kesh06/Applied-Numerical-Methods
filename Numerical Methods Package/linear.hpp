#ifndef NUMERICAL_METHODS_LINEAR_HPP
#define NUMERICAL_METHODS_LINEAR_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace NumericalMethods {
    namespace Linear {

        // --- Matrix Norms ---
        
        inline double norm_inf(const std::vector<std::vector<double>>& A) {
            double max_sum = 0.0;
            for (const auto& row : A) {
                double row_sum = 0.0;
                for (double val : row) row_sum += std::abs(val);
                if (row_sum > max_sum) max_sum = row_sum;
            }
            return max_sum;
        }

        inline double norm_l1(const std::vector<std::vector<double>>& A) {
            if (A.empty()) return 0.0;
            int n = A.size();
            int m = A[0].size();
            double max_sum = 0.0;
            for (int j = 0; j < m; j++) {
                double col_sum = 0.0;
                for (int i = 0; i < n; i++) col_sum += std::abs(A[i][j]);
                if (col_sum > max_sum) max_sum = col_sum;
            }
            return max_sum;
        }

        inline double norm_l2(const std::vector<std::vector<double>>& A, double eps = 1e-9, int max_iter = 1000) {
            int n = A.size();
            int m = A[0].size();
            std::vector<std::vector<double>> ATA(m, std::vector<double>(m, 0.0));
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < m; j++) {
                    for (int k = 0; k < n; k++) ATA[i][j] += A[k][i] * A[k][j];
                }
            }
            std::vector<double> v(m, 1.0);
            double lambda_old = 0.0, lambda_new = 0.0;
            for (int iter = 0; iter < max_iter; iter++) {
                std::vector<double> w(m, 0.0);
                for (int i = 0; i < m; i++) {
                    for (int j = 0; j < m; j++) w[i] += ATA[i][j] * v[j];
                }
                double numerator = 0.0, denominator = 0.0;
                for (int i = 0; i < m; i++) {
                    numerator += v[i] * w[i];
                    denominator += v[i] * v[i];
                }
                if (std::abs(denominator) < 1e-12) break;
                lambda_new = numerator / denominator;
                if (std::abs(lambda_new - lambda_old) < eps) break;
                lambda_old = lambda_new;
                double norm_w = 0.0;
                for (double val : w) norm_w += val * val;
                norm_w = std::sqrt(norm_w);
                if (norm_w < 1e-12) break;
                for (int i = 0; i < m; i++) v[i] = w[i] / norm_w;
            }
            return std::sqrt(std::max(0.0, lambda_new));
        }

        // --- Direct Solvers & Engine ---

        inline void to_REF(std::vector<std::vector<double>>& A) {
            if (A.empty()) return;
            int n = A.size(), m = A[0].size(), r = 0;
            for (int k = 0; k < m && r < n; k++) {
                int p = r;
                double max_val = std::abs(A[r][k]);
                for (int i = r + 1; i < n; i++) {
                    if (std::abs(A[i][k]) > max_val) {
                        max_val = std::abs(A[i][k]);
                        p = i;
                    }
                }
                if (std::abs(A[p][k]) < 1e-12) continue;
                if (p != r) std::swap(A[r], A[p]);
                for (int i = r + 1; i < n; i++) {
                    double factor = A[i][k] / A[r][k];
                    for (int j = k; j < m; j++) A[i][j] -= factor * A[r][j];
                }
                r++;
            }
        }

        inline std::vector<double> gauss(const std::vector<std::vector<double>>& A_, const std::vector<double>& b) {
            int n = A_.size();
            if (n != (int)A_[0].size()) throw std::runtime_error("Matrix must be square.");
            std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) A[i][j] = A_[i][j];
                A[i][n] = b[i];
            }
            to_REF(A);
            for (int i = 0; i < n; i++) {
                if (std::abs(A[i][i]) < 1e-12) throw std::runtime_error("Matrix is singular.");
            }
            std::vector<double> x(n, 0.0);
            x[n - 1] = A[n - 1][n] / A[n - 1][n - 1];
            for (int i = n - 2; i >= 0; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++) sum += A[i][j] * x[j];
                x[i] = (A[i][n] - sum) / A[i][i];
            }
            return x;
        }

        inline std::vector<double> gauss_scaled(const std::vector<std::vector<double>>& A_, const std::vector<double>& b) {
            int n = A_.size();
            if (n != (int)A_[0].size()) throw std::runtime_error("Matrix must be square.");
            std::vector<std::vector<double>> A(n, std::vector<double>(n + 1));
            std::vector<double> s(n, 0.0);
            for (int i = 0; i < n; i++) {
                double max_row_val = 0.0;
                for (int j = 0; j < n; j++) {
                    A[i][j] = A_[i][j];
                    if (std::abs(A[i][j]) > max_row_val) max_row_val = std::abs(A[i][j]);
                }
                A[i][n] = b[i];
                s[i] = max_row_val;
                if (std::abs(s[i]) < 1e-12) throw std::runtime_error("Matrix contains an all-zero row.");
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= n; j++) A[i][j] /= s[i];
            }
            to_REF(A);
            for (int i = 0; i < n; i++) {
                if (std::abs(A[i][i]) < 1e-12) throw std::runtime_error("Matrix is singular.");
            }
            std::vector<double> x(n, 0.0);
            x[n - 1] = A[n - 1][n] / A[n - 1][n - 1];
            for (int i = n - 2; i >= 0; i--) {
                double sum = 0.0;
                for (int j = i + 1; j < n; j++) sum += A[i][j] * x[j];
                x[i] = (A[i][n] - sum) / A[i][i];
            }
            return x;
        }

        // --- Iterative Core & Solvers ---

        inline std::vector<double> iterate_step(const std::vector<std::vector<double>>& M, const std::vector<double>& y_n, const std::vector<double>& c) {
            int n = M.size();
            std::vector<double> y_next(n, 0.0);
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < n; j++) sum += M[i][j] * y_n[j];
                y_next[i] = sum + c[i];
            }
            return y_next;
        }

        inline double compute_error(const std::vector<double>& x_next, const std::vector<double>& x) {
            double max_err = 0.0;
            for (size_t i = 0; i < x.size(); i++) {
                double err = std::abs(x_next[i] - x[i]);
                if (err > max_err) max_err = err;
            }
            return max_err;
        }

        inline std::vector<double> jacobi(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double eps, int max_iter = 1000) {
            int n = A.size();
            std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
            std::vector<double> c(n, 0.0);
            for (int i = 0; i < n; i++) {
                if (std::abs(A[i][i]) < 1e-12) throw std::runtime_error("Zero diagonal detected.");
                c[i] = b[i] / A[i][i];
                for (int j = 0; j < n; j++) {
                    if (i != j) M[i][j] = -A[i][j] / A[i][i];
                }
            }
            std::vector<double> x(n, 0.0);
            for (int iter = 0; iter < max_iter; iter++) {
                std::vector<double> x_next = iterate_step(M, x, c);
                if (compute_error(x_next, x) < eps) return x_next;
                x = x_next;
            }
            return x;
        }

        inline std::vector<double> gauss_seidel(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double eps, int max_iter = 1000) {
            int n = A.size();
            std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
            std::vector<double> c(n, 0.0);
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < i; j++) sum += A[i][j] * c[j];
                c[i] = (b[i] - sum) / A[i][i];
            }
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    double target = (i < j) ? -A[i][j] : 0.0;
                    double sum = 0.0;
                    for (int k = 0; k < i; k++) sum += A[i][k] * M[k][j];
                    M[i][j] = (target - sum) / A[i][i];
                }
            }
            std::vector<double> x(n, 0.0);
            for (int iter = 0; iter < max_iter; iter++) {
                std::vector<double> x_next = iterate_step(M, x, c);
                if (compute_error(x_next, x) < eps) return x_next;
                x = x_next;
            }
            return x;
        }

        inline std::vector<double> sor(const std::vector<std::vector<double>>& A, const std::vector<double>& b, double w, double eps, int max_iter = 1000) {
            int n = A.size();
            std::vector<std::vector<double>> M(n, std::vector<double>(n, 0.0));
            std::vector<double> c(n, 0.0);
            for (int i = 0; i < n; i++) {
                double sum = 0.0;
                for (int j = 0; j < i; j++) sum += w * A[i][j] * c[j];
                c[i] = (w * b[i] - sum) / A[i][i];
            }
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    double target = 0.0;
                    if (i == j)     target = (1.0 - w) * A[i][i];
                    else if (i < j) target = -w * A[i][j];
                    double sum = 0.0;
                    for (int k = 0; k < i; k++) sum += w * A[i][k] * M[k][j];
                    M[i][j] = (target - sum) / A[i][i];
                }
            }
            std::vector<double> x(n, 0.0);
            for (int iter = 0; iter < max_iter; iter++) {
                std::vector<double> x_next = iterate_step(M, x, c);
                if (compute_error(x_next, x) < eps) return x_next;
                x = x_next;
            }
            return x;
        }
    }
}

#endif

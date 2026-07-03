#ifndef NUMERICAL_METHODS_QUADRATURE_HPP
#define NUMERICAL_METHODS_QUADRATURE_HPP

#include <functional>
#include <stdexcept>

namespace NumericalMethods {
    namespace Quadrature {

        // =====================================================================
        // 1. Composite Trapezoidal Rule
        // =====================================================================
        inline double trapezoidal(std::function<double(double)> f, double a, double b, int n) {
            if (n <= 0) {
                throw std::invalid_argument("Error: Number of intervals (n) must be greater than 0.");
            }

            double h = (b - a) / n;
            double sum = 0.5 * (f(a) + f(b));

            for (int i = 1; i < n; i++) {
                double x = a + i * h;
                sum += f(x);
            }

            return sum * h;
        }

        // =====================================================================
        // 2. Composite Simpson's 1/3rd Rule
        // =====================================================================
        inline double simpson(std::function<double(double)> f, double a, double b, int n) {
            if (n <= 0) {
                throw std::invalid_argument("Error: Number of intervals (n) must be greater than 0.");
            }
            if (n % 2 != 0) {
                throw std::invalid_argument("Error: Number of intervals (n) must be an even integer for Simpson's rule.");
            }

            double h = (b - a) / n;
            double sum = f(a) + f(b);

            for (int i = 1; i < n; i++) {
                double x = a + i * h;
                if (i % 2 == 0) {
                    sum += 2.0 * f(x); // Even indices get a weight of 2
                } else {
                    sum += 4.0 * f(x); // Odd indices get a weight of 4
                }
            }

            return (h / 3.0) * sum;
        }
    }
}

#endif
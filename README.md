# Numerical Methods Suite

This repository contains numerical analysis implementations developed for university coursework, featuring both individual raw algorithmic scripts and a production-grade C++ library.

---

## Repository Structure

* **`Lab Assignments/`**: Standalone MATLAB scripts developed for weekly laboratory assignments, implementing core numerical routines.
* **`Lab Exam/`**: Solution set and implementation scripts used during the final laboratory examination.
* **`Numerical Methods Package/`**: A modular, STL-style C++ library providing optimized, production-ready headers for mathematical computing.

---

## C++ Library Modules

The custom C++ implementation is organized into stateless, decoupled namespaces:

* **`Linear`**: Features automated row scaling, partial pivoting Gaussian elimination, matrix norm engines ($L_1$, $L_{\infty}$, and iterative $L_2$ via Rayleigh Quotient), and iterative solvers (Jacobi, Gauss-Seidel, SOR).
* **`Interpolation`**: Includes Newton’s Divided Difference, Lagrange Polynomials, and Natural Cubic Splines integrated with the internal linear system engine.
* **`NonLinear`**: Implements robust root-finding routines including Bisection, Regula Falsi, Secant, and Newton-Raphson methods with low-derivative protections.
* **`Quadrature`**: Provides definitive numeric integration interfaces implementing the Composite Trapezoidal and Composite Simpson’s 1/3rd rules.

---

## Example C++ Library Usage

Compile `main.cpp` using any standard C++17 compiler to execute the comprehensive verification benchmark suite:

```bash
g++ -std=c++17 main.cpp -o numerical_bench
./numerical_bench

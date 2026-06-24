# Numerical Analysis

This repository contains MATLAB implementations of key numerical analysis algorithms. It includes methods for finding roots of nonlinear equations, as well as direct and iterative methods for solving systems of linear equations.

## File Description & Features

### 1. Root Finding for Nonlinear Equations
* **`Bisection_Newton-Raphson_compare.m`**
  * A script that compares the performance and convergence processes of the **Bisection Method** and the **Newton-Raphson Method**.
  * It demonstrates the differences between the stable Bisection method, which repeatedly halves intervals, and the fast-converging Newton-Raphson method, which utilizes derivatives.

### 2. Solving Systems of Linear Equations
* **`System_of_Linear_Equations.m`**
  * A comprehensive testing script for solving systems of linear equations.
* **`naive_gauss.m`**
  * Implementation of **Naive Gauss Elimination**. It solves the equations directly through forward elimination and backward substitution without partial pivoting.
* **`jacobi_method.m`**
  * Implementation of the **Jacobi Method**. Starting from an initial guess, it iteratively calculates the next approximation based on all values from the previous step.
* **`gauss_seidel.m`**
  * Implementation of the **Gauss-Seidel Method**. Unlike the Jacobi method, it uses the most recently calculated (updated) values immediately in the calculation of the next variable, which generally improves the convergence speed.

## Language & Environment
* **Language:** MATLAB

## How to Use
1. Clone this repository to your local machine:
```bash
   git clone [https://github.com/CNchanhoe/Numerical_analysis.git](https://github.com/CNchanhoe/Numerical_analysis.git)

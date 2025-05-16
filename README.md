# Numerical optimization for large scale problems
This repository contains the code for the test of convergence of three different functions using both Nelder Mead and Modified Newton methods.
Nelder Mead is performed with $n=10,25,50$ dimensions.
Modified Newton is performed with $n=1e^3$, $1e^4$, $1e^5$ dimensions.

## Nelder-Mead Method

- **Reflection**: $x_r = \bar{x} + \alpha(\bar{x} - x_h)$
- **Expansion**:  $x_e = \bar{x} + \gamma(x_r - \bar{x})$
- **Contraction**: $x_c = \bar{x} + \rho(x_h - \bar{x})$
- **Reduction**: $x_i = x_l + \sigma(x_i - x_l), \quad \forall i \neq l$

Where:
- $\bar{x}$ is the centroid of the simplex excluding the worst point $x_h$
- $\alpha > 0$, $\gamma > 1$, $0 < \rho < 0.5$, and $0 < \sigma < 1$ are algorithm parameters

---

## Modified Newton Method

The Modified Newton Method is a robust variant of the classical Newton's method for optimization. It uses second-order information (the Hessian matrix) but modifies it to ensure a **descent direction**, improving convergence when the Hessian is not positive definite.

### Update formula:

$$
x_{k+1} = x_k - \alpha_k H_k^{-1} \nabla f(x_k)
$$

Where:
- $\nabla f(x_k)$ is the gradient of the function
- $H_k$ is the Hessian matrix (possibly modified to be positive definite)
- $\alpha_k$ is the step size found through line search (e.g., Armijo rule)

The modification of the Hessian is typically done by adding a diagonal shift:

$$
H_k^{\text{mod}} = H_k + \mu I
$$

with $\mu > 0$ chosen such that $H_k^{\text{mod}} \succ 0$.

---

## How to Run

To run either method, navigate to the respective folder (`Nelder-Mead__method/` or `modified_newton_method/`) and execute the `main.py` script:

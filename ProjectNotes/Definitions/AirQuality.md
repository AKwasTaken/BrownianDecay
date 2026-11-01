# Derivation of the Itô SDE for Indoor Air Quality

## Step 1: The Deterministic Mass-Balance Equation
We begin with the established ordinary differential equation (Equation 6.9) governing the indoor concentration of the pollutant:
$$v \dot{C}^i = -k \left(q^0 + f^1 q^1 + q^2\right) C^i + S - R + k \left(q^0(1 - f^0) + q^2\right)C^0$$
## Step 2: Define the Stochastic Parameters
To inject real-world fluctuations, we decompose the four highly variable parameters into a deterministic mean (denoted with a bar) and a stochastic noise term (denoted by $\xi$):
$$S = \bar{S} + \xi^{(1)}$$
$$C^0 = \bar{C}^0 + \xi^{(2)}$$
$$q^0 = \bar{q}^0 + \xi^{(3)}$$
$$q^2 = \bar{q}^2 + \xi^{(4)}$$
## Step 3: Substitution
We substitute these decomposed parameters back into the baseline equation:
$$v \dot{C}^i = -k \left( (\bar{q}^0 + \xi^{(3)}) + f^1 q^1 + (\bar{q}^2 + \xi^{(4)}) \right) C^i \\ + (\bar{S} + \xi^{(1)}) - R \\ + k \left( (\bar{q}^0 + \xi^{(3)})(1 - f^0) + (\bar{q}^2 + \xi^{(4)}) \right) (\bar{C}^0 + \xi^{(2)})$$
## Step 4: Expansion and Distribution
Next, we expand the brackets to separate the deterministic terms from the noise terms.
**First term (Losses):**
$$-k(\bar{q}^0 + f^1 q^1 + \bar{q}^2)C^i - k(\xi^{(3)} + \xi^{(4)})C^i$$
**Second term (Sources & Sinks):**
$$\bar{S} - R + \xi^{(1)}$$
**Third term (Outdoor Contributions):** Expand the product of the flow rates and the outdoor concentration:
$$k \left[ \bar{q}^0(1 - f^0) + \xi^{(3)}(1 - f^0) + \bar{q}^2 + \xi^{(4)} \right] \left[ \bar{C}^0 + \xi^{(2)} \right]$$
Multiplying this out yields four distinct components:
1. $k (\bar{q}^0(1 - f^0) + \bar{q}^2) \bar{C}^0$ _(Purely deterministic)_
2. $k (\bar{q}^0(1 - f^0) + \bar{q}^2) \xi^{(2)}$ _(Noise from outdoor concentration)_
3. $k (\xi^{(3)}(1 - f^0) + \xi^{(4)}) \bar{C}^0$ _(Noise from flow rates)_
4. $k (\xi^{(3)}(1 - f^0) + \xi^{(4)}) \xi^{(2)}$ _(Cross-product of noise terms)_
## Step 5: Neglecting Higher-Order Noise
In stochastic calculus, the product of two independent white-noise fluctuations (like $\xi^{(2)}\xi^{(3)}$ or $\xi^{(2)}\xi^{(4)}$) is infinitesimally small and mathematically collapses when transitioning to the Itô integral limit. Therefore, we drop the cross-product term:
$$k (\xi^{(3)}(1 - f^0) + \xi^{(4)}) \xi^{(2)} \approx 0$$
## Step 6: Grouping Deterministic (Drift) and Stochastic (Diffusion) Terms
We now gather all the purely deterministic terms to form the drift function $g(C^i, t)$, and factor out the individual noise terms $\xi^{(1)}, \xi^{(2)}, \xi^{(3)}, \xi^{(4)}$ to form the diffusion vector $H(C^i, t)$.
**The Drift Term $g(C^i, t)$:**
$$g(C^i, t) = -k(\bar{q}^0 + f^1 q^1 + \bar{q}^2)C^i + \bar{S} - R + k(\bar{q}^0(1 - f^0) + \bar{q}^2)\bar{C}^0$$
**The Diffusion Term $H(C^i, t)$:** Grouping the coefficients for each $\xi$:
- For $\xi^{(1)}$: $1$
- For $\xi^{(2)}$: $k(\bar{q}^0(1 - f^0) + \bar{q}^2)$
- For $\xi^{(3)}$: $k\bar{C}^0(1 - f^0) - kC^i = k(\bar{C}^0(1 - f^0) - C^i)$
- For $\xi^{(4)}$: $k\bar{C}^0 - kC^i = k(\bar{C}^0 - C^i)$
> **Note on a likely typo in the source text:** The provided text lists the third row of the $H$ vector as $k(\bar{q}^0(1 - f^0) - C^i)\bar{C}^0$. Based on the step-by-step algebraic expansion above, the correct coefficient for $\xi^{(3)}$ is $k(\bar{C}^0(1 - f^0) - C^i)$.

Writing these coefficients as a column vector (the transpose of the row vector $H$):
$$H(C^i, t)^T = \begin{pmatrix} 1 \\ k \left(\bar{q}^0 (1 - f^0) + \bar{q}^2\right) \\ k(\bar{C}^0(1 - f^0) - C^i) \\ k(\bar{C}^0 - C^i) \end{pmatrix}$$
## Step 7: Final SDE Formulation
By multiplying the entire equation by $dt$ and converting the noise terms $\xi_t dt$ into the increments of independent Wiener processes $d\xi_t$, we arrive at the final Itô Stochastic Differential Equation (Equation 6.10):
$$v dC_t^i = g(C^i, t) dt + H(C^i, t) d\xi_t$$
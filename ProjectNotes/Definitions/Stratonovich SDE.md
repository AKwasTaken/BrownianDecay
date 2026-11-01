# The Stratonovich vs. Itô Interpretations in Stochastic Modelling
In stochastic calculus, integrating a differential equation driven by a Wiener process (Brownian motion), $W_t$, presents a unique mathematical ambiguity. Because continuous white noise fluctuates infinitely fast and has unbounded variation, the value of the stochastic integral depends entirely on where the state is evaluated within a given infinitesimal time step $[t, t+dt]$. This ambiguity yields two primary frameworks: the Itô and Stratonovich interpretations.
#### **1. Mathematical Interpretations and the $\circ$ Notation**
- **The Itô Interpretation ($dW_t$):** Evaluates the state variable strictly at the left endpoint (the beginning) of the time step. This means the system state is entirely uncorrelated with the noise increment that occurs over that step. Itô calculus preserves the non-anticipating (martingale) property, making it standard in financial mathematics, but it requires a modified chain rule (Itô's Lemma).
- **The Stratonovich Interpretation ($\circ dW_t$):** Evaluates the state variable at the midpoint of the time step. The $\circ$ symbol explicitly denotes this midpoint evaluation. Because it allows a slight correlation between the state and the noise increment, Stratonovich calculus obeys the standard Newton-Leibniz rules of classical calculus (the ordinary chain rule applies)
#### **2. The Conversion Formulation**
Because the Stratonovich midpoint evaluation implicitly correlates the state with the noise, it generates a systematic bias over time. When converting to Itô rules (which assume zero correlation), this bias must be explicitly added to the deterministic drift of the equation. This correction is known as the "noise-induced drift."
For a general Stratonovich SDE with a deterministic rate $F(X_t)$ and multiple independent noise sources $g_i(X_t)$:

$$dX_t = F(X_t)dt + \sum_{i} g_i(X_t) \circ dW_t^i$$
The equivalent Itô SDE is obtained by adding the derivative-based correction term $\frac{1}{2} g_i(X_t) g_i'(X_t)$ for each multiplicative noise term:
$$dX_t = \left[ F(X_t) + \frac{1}{2} \sum_{i} g_i(X_t) \frac{\partial g_i(X_t)}{\partial X_t} \right] dt + \sum_{i} g_i(X_t) dW_t^i$$
By applying this transformation, the $\circ$ symbol is dropped, indicating that standard Itô numerical integration can now be safely applied without violating the underlying physics of the system.

| **Feature**          | **Itô SDE**                                                                          | **Stratonovich SDE (∘)**                                                                                    |
| -------------------- | ------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------- |
| **Notation**         | $g(X_t) \, dW_t$ (no symbol)                                                         | $g(X_t) \circ dW_t$ (uses $\circ$)                                                                          |
| **Evaluation Point** | Left endpoint ($t_i$)                                                                | Midpoint ($\frac{t_i + t_{i+1}}{2}$)                                                                        |
| **Calculus Rules**   | Requires **Itô's Lemma** (adds a second-derivative correction term)                  | Follows **Ordinary Calculus** (standard chain rule works)                                                   |
| **Physical Meaning** | Preferred in **Finance/Economics** (where future info at $t_{i+1}$ isn't available). | Preferred in **Physics/Engineering** (when white noise is an approximation of real-world continuous noise). |

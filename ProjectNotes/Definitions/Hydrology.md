# Derivation for Stratonovich SDE
### **Step 1: Write out the deterministic rate equation**
Start with the rate of change of the reservoir outflow over time $\frac{dX_t}{dt} = F(X_t)$ using Equation (6.7):

$$
\frac{dX_t}{dt} = AP X_t^\alpha - AP\Phi(1 - \psi(X_t)) X_t^\alpha - AE\psi(X_t)X_t^\alpha - AX_t^\beta - BX_t^\gamma
$$

### **Step 2: Substitute the noisy environmental variables**
Replace the environmental terms $P$, $E$, and $\Phi$ with their noisy representations:
- $P \rightarrow P + \sigma_1 \xi_t^1$
- $E \rightarrow E + \sigma_2 \xi_t^2$
- $\Phi \rightarrow \Phi + \sigma_3 \xi_t^3$
Substituting these into the rate equation gives:

$$
\begin{aligned} \frac{dX_t}{dt} = & A(P + \sigma_1 \xi_t^1) X_t^\alpha \\ & - A(P + \sigma_1 \xi_t^1)(\Phi + \sigma_3 \xi_t^3)(1 - \psi(X_t)) X_t^\alpha \\ & - A(E + \sigma_2 \xi_t^2)\psi(X_t)X_t^\alpha \\ & - AX_t^\beta - BX_t^\gamma \end{aligned}
$$

### **Step 3: Expand the product in the second term**
Expanding the second term $(P + \sigma_1 \xi_t^1)(\Phi + \sigma_3 \xi_t^3)$ yields four components:

$$
(P + \sigma_1 \xi_t^1)(\Phi + \sigma_3 \xi_t^3) = P\Phi + P\sigma_3\xi_t^3 + \Phi\sigma_1\xi_t^1 + \sigma_1\sigma_3\xi_t^1\xi_t^2
$$

Because $\xi_t^1$ and $\xi_t^3$ are independent white noise processes, their product term $\sigma_1\sigma_3\xi_t^1\xi_t^3$ has a mean of zero and negligible variance relative to first-order noise terms, so it is dropped. This simplifies the expansion to:

$$
P\Phi + P\sigma_3\xi_t^3 + \Phi\sigma_1\xi_t^1
$$

### **Step 4: Expand all terms fully**
Substitute the simplified product back into the full differential equation:

$$
\begin{aligned} \frac{dX_t}{dt} = & AP X_t^\alpha + A\sigma_1\xi_t^1 X_t^\alpha \\ & - AP\Phi(1 - \psi(X_t))X_t^\alpha - A P \sigma_3 \xi_t^3 (1 - \psi(X_t))X_t^\alpha - A \Phi \sigma_1 \xi_t^1 (1 - \psi(X_t))X_t^\alpha \\ & - AE\psi(X_t)X_t^\alpha - A\sigma_2\xi_t^2\psi(X_t)X_t^\alpha \\ & - AX_t^\beta - BX_t^\gamma \end{aligned}
$$

### **Step 5: Separate and group deterministic vs. stochastic terms**
#### 1. **Deterministic Group (baseline $F(X_t)$):**
Collect all terms that do not contain $\xi_t$:

$$
AP X_t^\alpha - AP\Phi(1 - \psi(X_t))X_t^\alpha - AE\psi(X_t)X_t^\alpha - AX_t^\beta - BX_t^\gamma = F(X_t)
$$

#### 2. **Noise 1 Group ($\xi_t^1$ terms):**
Collect all terms containing $\sigma_1 \xi_t^1$:

$$
A\sigma_1\xi_t^1 X_t^\alpha - A \Phi \sigma_1 \xi_t^1 (1 - \psi(X_t))X_t^\alpha
$$

Factor out $\sigma_1 \xi_t^1$:

$$
= \sigma_1 \left( A X_t^\alpha - A\Phi(1 - \psi(X_t))X_t^\alpha \right) \xi_t^1
$$

Factor out $A X_t^\alpha$:

$$
= \sigma_1 \left( 1 - \Phi(1 - \psi(X_t)) \right) A X_t^\alpha \xi_t^1
$$

Expand inside the brackets:

$$
= \sigma_1 \left( 1 - \Phi + \Phi\psi(X_t) \right) A X_t^\alpha \xi_t^1 \quad \text{or} \quad \sigma_1 \left( 1 - \Phi + A\Phi\psi(X_t) X_t^\alpha \right) \xi_t^1
$$

#### 3. **Noise 2 Group ($\xi_t^2$ terms):**
Collect terms containing $\sigma_2 \xi_t^2$:

$$
- A\sigma_2\xi_t^2\psi(X_t)X_t^\alpha = - \sigma_2 A \psi(X_t) \xi_t^2
$$

#### 4. **Noise 3 Group ($\xi_t^3$ terms):**
Collect terms containing $\sigma_3 \xi_t^3$:

$$
- A P \sigma_3 \xi_t^3 (1 - \psi(X_t))X_t^\alpha = - \sigma_3 A P (1 - \psi(X_t)) X_t^\alpha \xi_t^3
$$

Combining these groups gives **Equation 1**:

$$
\frac{dX_t}{dt} = F(X_t) + \sigma_1 \left( 1 - \Phi + A\Phi\psi(X_t)X_t^\alpha \right) \xi_t^1 - \sigma_2 A \psi(X_t) \xi_t^2 - \sigma_3 A P (1 - \psi(X_t)) X_t^\alpha \xi_t^3
$$

### **Step 6: Convert from differential rate form to SDE form (Equation 2)**
1. **Multiply the entire equation by $dt$:**
    
$$
\left(\frac{dX_t}{dt}\right) dt = F(X_t)dt + g_1(X_t)\left(\xi_t^1 dt\right) - g_2(X_t)\left(\xi_t^2 dt\right) - g_3(X_t)\left(\xi_t^3 dt\right)
$$

2. **Apply Langevin's definition of Gaussian white noise:** In stochastic calculus, continuous physical white noise $\xi_t$ is defined as the time derivative of a standard Wiener process (Brownian motion) $W_t$:
    
$$
\xi_t = \frac{dW_t}{dt} \implies \xi_t dt = dW_t
$$

3. **Adopt Stratonovich notation ($\circ dW_t$):** Because physical white noise fluctuations in environmental systems (like rainfall and evaporation) are continuous approximations of real-world physical processes, the integral obeys standard classical calculus rules. This midpoint evaluation rule is denoted by the Stratonovich symbol $\circ dW_t$.

Replacing $\left(\frac{dX_t}{dt}\right)dt$ with $dX_t$ and each $\xi_t^i dt$ with $\circ dW_t^i$ yields **Equation 2**:

$$
dX_t = F(X_t)dt + \sigma_1(1 - \Phi + A\Phi\psi(X_t)X_t^\alpha) \circ dW_t^1 - \sigma_2 A \psi(X_t) \circ dW_t^2 - \sigma_3 A P (1 - \psi(X_t))X_t^\alpha \circ dW_t^3
$$

____
# Converting the Stratonovich SDE to an Itô SDE
If we want to simulate this system using standard numerical methods (like the Euler-Maruyama method) or apply Itô calculus, we must convert the Stratonovich SDE into an Itô SDE. We cannot simply drop the Stratonovich circle ($\circ$) because Stratonovich integration evaluates the state at the midpoint of a time step (preserving classical chain-rule calculus for physical processes), while Itô integration evaluates the state at the beginning of the time step.
To switch to Itô notation, we must add a "noise-induced drift" correction term to account for how the multiplicative noise interacts with the state variable over the time step.
#### **Step A: Define the corrected diffusion terms**
First, we correct a minor algebraic reduction from Step 5 to ensure the state variable $X_t^\alpha$ is properly retained across all noise terms. We define our three noise coefficients $g_1(X_t)$, $g_2(X_t)$, and $g_3(X_t)$ as:
- **Noise 1 ($\xi_t^1$):** $g_1(X_t) = A \sigma_1 X_t^\alpha \big(1 - \Phi + \Phi\psi(X_t)\big)$
- **Noise 2 ($\xi_t^2$):** $g_2(X_t) = -A \sigma_2 X_t^\alpha \psi(X_t)$
- **Noise 3 ($\xi_t^3$):** $g_3(X_t) = -A \sigma_3 P X_t^\alpha \big(1 - \psi(X_t)\big)$
#### **Step B: Calculate the Itô drift correction**
The general formula to convert a Stratonovich SDE $dX_t = F(X_t)dt + \sum g_i(X_t) \circ dW_t^i$ to an Itô SDE is:

$$
dX_t = \left( F(X_t) + \frac{1}{2}\sum_{i} g_i(X_t) \frac{\partial g_i(X_t)}{\partial X_t} \right) dt + \sum_{i} g_i(X_t) dW_t^i
$$

To find the correction factor $\frac{1}{2} g_i(X_t) g_i'(X_t)$ for each term, we take the derivative of each diffusion term with respect to $X_t$, applying the product rule since both $X_t^\alpha$ and $\psi(X_t)$ depend on $X_t$:
1. **Derivative of $g_1(X_t)$:**

$$
g_1'(X_t) = A \sigma_1 \left[ \alpha X_t^{\alpha-1} \big(1 - \Phi + \Phi\psi(X_t)\big) + X_t^\alpha \Phi \psi'(X_t) \right]
$$

2. **Derivative of $g_2(X_t)$:**

$$
g_2'(X_t) = -A \sigma_2 \left[ \alpha X_t^{\alpha-1} \psi(X_t) + X_t^\alpha \psi'(X_t) \right]
$$

3. **Derivative of $g_3(X_t)$:**

$$
g_3'(X_t) = -A \sigma_3 P \left[ \alpha X_t^{\alpha-1} \big(1 - \psi(X_t)\big) - X_t^\alpha \psi'(X_t) \right]
$$

#### **Step C: Write the final Itô SDE**
By adding the combined drift correction to the original deterministic rate equation $F(X_t)$, we arrive at the standard Itô SDE (note the absence of the $\circ$ symbol):

$$
dX_t = \left[ F(X_t) + \frac{1}{2} \left( g_1(X_t)g_1'(X_t) + g_2(X_t)g_2'(X_t) + g_3(X_t)g_3'(X_t) \right) \right] dt + g_1(X_t) dW_t^1 + g_2(X_t) dW_t^2 + g_3(X_t) dW_t^3
$$


_(To write this equation out fully in a paper, you would substitute the $g_i$ and $g_i'$ expressions defined in Steps A and B into this final structure)._
____

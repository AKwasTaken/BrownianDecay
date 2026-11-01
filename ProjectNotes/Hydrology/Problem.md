# The Core Problem & The Physical Setup
The goal here is to predict the **outflow** of water from a single lake or reservoir.
If you are managing a dam or a reservoir, you need to know exactly how much water is going to flow out of it. However, the amount of water in a lake isn't just about what you pour in; it's a massive natural balancing act influenced by:
- **Precipitation ($P$):** Rain falling on the lake and the surrounding land (catchment area).
- **Evaporation ($E$):** Water evaporating off the surface into the air.
- **Seepage:** Water quietly leaking out through the soil/ground at the bottom.
- **Abstractions ($\Phi$):** Other environmental factors or water being drawn out.
____
# The Baseline Math (The Perfect World)
Before adding randomness, hydrologists used a standard **Ordinary Differential Equation (ODE)**(similar to Population modelling) to model this:
$$(6.6) \qquad \dot{x} = F(x)$$
This simply says: "The rate at which the outflow changes ($\dot{x}$) is a function of current conditions ($F(x)$)."
Equation 6.7 defines exactly what that function is:
$$(6.7) \qquad F(x) = AP x^\alpha - AP\Phi(1 - \psi(x)) x^\alpha - AE\psi(x)x^\alpha - Ax^\beta - Bx^\gamma$$
**Breaking down the alphabet soup of 6.7:**
- **$x$**: The outflow of the lake.
- **The $P$ terms**: Represent how rainfall adds to the water.
- **The $E$ term**: Represents how evaporation subtracts water.
- **$A$:** A geographic constant representing the surface area of the lake and its surrounding catchment basin.
- **$\Phi$ (Phi - Abstractions):** This is usually a fractional constant (a number between 0 and 1) representing the percentage of water artificially pumped out for human use, like agriculture, industrial cooling, or municipal drinking water.
- **$Bx^\gamma$ (the last term)**: Represents the seepage—water escaping into the ground.
- **$\psi(x)$**: This represents the ratio of the surface area of the surrounding land (catchment) to the surface area of the reservoir itself.
- **$A, B, C, D, \alpha, \beta, \gamma, \delta$**: These are just structural constants and exponents used to fit the equation to a specific real-world lake's shape and soil type.
Ultimately, Equation 6.7 is just a highly detailed, localized way of mathematically writing:
**F(x) = (Rainfall entering) - (Water pumped out by humans) - (Water evaporating) - (Main dam outflow) - (Water leaking into the ground).**

| **Mathematical Term in Equation 6.7** | **Real-World Physical Meaning**                                                             | **Key Variables & Constants**                                                                  |
| ------------------------------------- | ------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------- |
| **$AP x^\alpha$**                     | **Rainfall entering** (Gross precipitation falling on the lake and surrounding land)        | **$A$**: Catchment area<br>**$P$**: Precipitation rate<br>**$\alpha$**: Shape scaling exponent |
| **$- AP\Phi(1 - \psi(x)) x^\alpha$**  | **Water pumped out by humans** (Abstractions for agriculture, industry, or municipal use)   | **$\Phi$**: Percentage pumped out<br>**$\psi(x)$**: Surface area ratio                         |
| **$- AE\psi(x)x^\alpha$**             | **Water evaporating** (Water vaporizing off the lake surface into the atmosphere)           | **$E$**: Evaporation rate<br>**$\psi(x)$**: Surface area ratio                                 |
| **$- Ax^\beta$**                      | **Main dam outflow** (The primary water released down the river or over the spillway)       | **$\beta$**: Shape scaling exponent for main outflow                                           |
| **$- Bx^\gamma$**                     | **Water leaking into the ground** (Seepage lost through porous soil and rock at the bottom) | **$B$**: Seepage coefficient<br>**$\gamma$**: Shape scaling exponent for depth/pressure        |
![[Pasted image 20261101190849.png]]
____
# Injecting the Real-World Noise
The equation above assumes we know exactly how much it will rain ($P$), evaporate ($E$), or be absorbed ($\Phi$). But weather is inherently unpredictable.
Researchers **Unny and Karmeshu** argued that these shouldn't be fixed numbers. Instead, they should be treated as **random variables**—an average expected value, plus a random, unpredictable fluctuation (Gaussian white noise).
They replaced the static variables with noisy ones:
- **Rainfall ($P$)** becomes: Average Rain ($P$) + Random Rain Fluctuation ($\sigma_1 \xi_t^1$)
- **Evaporation ($E$)** becomes: Average Evaporation ($E$) + Random Evap Fluctuation ($\sigma_2 \xi_t^2$)
- **Abstractions ($\Phi$)** becomes: Average Abstractions ($\Phi$) + Random Fluctuation ($\sigma_3 \xi_t^3$)
_(Note: The $\sigma$ values determine how wild the weather swings are, and $\xi_t$ represents the actual random "static" or noise)._
____
# The Final Stochastic Model (Equation 6.8)
They took these new, noisy variables and shoved them back into the original equation. _(The text mentions substituting this into "4.2", which is almost certainly a typo by the author and is meant to refer to 6.7)._
When you do the calculus to separate the predictable stuff from the random stuff, you get Equation 6.8:

$$(6.8) \qquad dX_t = F(X_t)dt + \sigma_1(1 - \Phi + A\Phi\psi(X_t)(X_t)^\alpha) \circ dW_t^1 \\ \qquad \qquad - \sigma_2 A \psi(X_t) \circ dW_t^2 \\ \qquad \qquad - \sigma_3 A P (1 - \psi(X_t))(X_t)^\alpha \circ dW_t^3$$

**What this means:**
- **$dX_t$**: The change in the lake's outflow over time.
- **$F(X_t)dt$**: The deterministic, baseline behaviour (this is just equation 6.7 acting normally).
- **The $dW_t$ terms**: These are the "Wiener processes" (standard mathematical models for random motion/noise).
    - $dW_t^1$ is the rainfall noise.
    - $dW_t^2$ is the evaporation noise (notice it has a minus sign, because evaporation removes water).
    - $dW_t^3$ is the abstraction noise (also subtracted).
**A quick math note:** You might notice the little circle ($\circ$) before the $dW_t$ terms. In the previous wastewater model, they used "Itô" calculus. Here, they use **Stratonovich** calculus. Both are just different mathematical rulebooks for how to integrate random noise. Physicists and engineers often prefer Stratonovich for physical systems driven by real-world "coloured" noise (like weather).
____
## Formulation of the Stochastic Differential Equation  

### Step 1: The Base Differential Equation
In calculus, the rate of change of water volume over time $t$ is written as $\frac{dX_t}{dt}$.
From equation (6.6):
$$\frac{dX_t}{dt} = F(X_t)$$
If you write out $F(X_t)$ fully using equation (6.7), it is:
$$\frac{dX_t}{dt} = A P X_t^\alpha - A P \Phi (1 - \psi(X_t)) X_t^\alpha - A E \psi(X_t) X_t^\alpha - A X_t^\beta - B X_t^\gamma$$
### Step 2: Substituting the Noisy Parameters
Unny and Karmeshu model the three environmental parameters as baseline means plus Gaussian white noise ($\xi_t$):
- $P \rightarrow P + \sigma_1 \xi_t^1$
- $E \rightarrow E + \sigma_2 \xi_t^2$
- $\Phi \rightarrow \Phi + \sigma_3 \xi_t^3$
Plugging these exact replacements into the $\frac{dX_t}{dt}$ equation gives:
$$\frac{dX_t}{dt} = A(P + \sigma_1 \xi_t^1)X_t^\alpha - A(P + \sigma_1 \xi_t^1)(\Phi + \sigma_3 \xi_t^3)(1 - \psi(X_t))X_t^\alpha - A(E + \sigma_2 \xi_t^2)\psi(X_t)X_t^\alpha - AX_t^\beta - BX_t^\gamma$$
### Step 3: Expanding and Grouping the Terms
If you expand those brackets algebraically, you get two distinct categories of terms:
1. **Deterministic terms:** The terms containing just $P$, $E$, and $\Phi$ (with no $\xi_t$). These collect together to form the original baseline function $F(X_t)$.
2. **Stochastic terms:** The terms containing the noise signals $\xi_t^1$, $\xi_t^2$, and $\xi_t^3$.
_(Note on the product term: As mentioned in your excerpt, the product of two noisy inputs $(P + \sigma_1\xi_t^1)(\Phi + \sigma_3\xi_t^3)$ gives a cross-term $\sigma_1\sigma_3 \xi_t^1 \xi_t^3$. Because $\xi_t^1$ and $\xi_t^3$ are independent white noise processes, their product is mathematically negligible and dropped, leaving $P\Phi + P\sigma_3\xi_t^3 + \Phi\sigma_1\xi_t^1$.)_
After grouping, the rate equation looks like this:
$$\frac{dX_t}{dt} = F(X_t) + \sigma_1 \left(1 - \Phi + A\Phi\psi(X_t)X_t^\alpha\right)\xi_t^1 - \sigma_2 A \psi(X_t)\xi_t^2 - \sigma_3 A P (1 - \psi(X_t))X_t^\alpha \xi_t^3$$
### Step 4: The Stochastic Calculus Step (Converting $\xi_t dt$ to $dW_t$)
This is where standard differential equations transform into Stochastic Differential Equations.
In stochastic modelling, formal "Gaussian white noise" ($\xi_t$) is defined as the time-derivative of a Brownian motion/Wiener process ($W_t$):
$$\xi_t = \frac{dW_t}{dt}$$
If you multiply the entire equation through by $dt$:
$$\left(\frac{dX_t}{dt}\right) dt = F(X_t)dt + g_1(X_t)\left(\xi_t^1 dt\right) - g_2(X_t)\left(\xi_t^2 dt\right) - g_3(X_t)\left(\xi_t^3 dt\right)$$
Since $\xi_t dt = \left(\frac{dW_t}{dt}\right) dt = dW_t$:
1. $\frac{dX_t}{dt} \cdot dt$ becomes **$dX_t$**
2. $F(X_t) \cdot dt$ becomes **$F(X_t)dt$**
3. Each noise term $\xi_t^i dt$ becomes **$dW_t^i$**
Because physical system parameters like rain and evaporation have non-zero correlation times in real life (coloured noise limit), this physical white noise approximation naturally yields a **Stratonovich SDE**, marked by the circle symbol ($\circ dW_t$):
$$dX_t = F(X_t)dt + \sigma_1(1 - \Phi + A\Phi\psi(X_t)X_t^\alpha) \circ dW_t^1 - \sigma_2 A \psi(X_t) \circ dW_t^2 - \sigma_3 A P (1 - \psi(X_t))X_t^\alpha \circ dW_t^3$$

**More in-depth derivation**:- [[Hydrology#Derivation for Stratonovich SDE|Derivation]]
___
# Modelling
While the Stratonovich interpretation correctly preserves classical calculus rules for physical environmental variables, directly simulating Stratonovich SDEs is computationally challenging. To utilize standard numerical solvers like the Euler-Maruyama method, we convert the system to its equivalent Itô form by introducing a noise-induced drift correction.
The step-by-step mathematical conversion is detailed here: [[Hydrology#Converting the Stratonovich SDE to an Itô SDE|Derivation]].
Applying this transformation yields the final Itô SDE used for our numerical modelling:
$$dX_t = \left[ F(X_t) + \frac{1}{2} \sum_{i=1}^{3} g_i(X_t)g_i'(X_t) \right] dt + \sum_{i=1}^{3} g_i(X_t) dW_t^i$$
where the multiplicative diffusion terms $g_i(X_t)$ for each noise source are defined as:
$$g_1(X_t) = A \sigma_1 X_t^\alpha \big(1 - \Phi + \Phi\psi(X_t)\big)$$
$$g_2(X_t) = -A \sigma_2 X_t^\alpha \psi(X_t)$$
$$g_3(X_t) = -A \sigma_3 P X_t^\alpha \big(1 - \psi(X_t)\big)$$
and their corresponding derivatives $g_i'(X_t) = \frac{\partial g_i}{\partial X_t}$, required to calculate the drift correction, are:
$$g_1'(X_t) = A \sigma_1 \left[ \alpha X_t^{\alpha-1} \big(1 - \Phi + \Phi\psi(X_t)\big) + X_t^\alpha \Phi \psi'(X_t) \right]$$
$$g_2'(X_t) = -A \sigma_2 \left[ \alpha X_t^{\alpha-1} \psi(X_t) + X_t^\alpha \psi'(X_t) \right]$$
$$g_3'(X_t) = -A \sigma_3 P \left[ \alpha X_t^{\alpha-1} \big(1 - \psi(X_t)\big) - X_t^\alpha \psi'(X_t) \right]$$
_______
# The Bigger Picture: The Great Lakes
Unny and Karmeshu used this equation in computer simulations to generate "sample paths"—essentially running hundreds of virtual simulations of a lake's outflow to see all the different ways the water levels might behave under random weather conditions.
But they didn't stop at one lake. They realized you could chain these equations together to model a **cascade of lakes**, like the Great Lakes in North America. In a cascade, the unpredictable outflow ($X_t$) of Lake Superior becomes a highly chaotic input for Lake Huron, which adds its own random weather noise, and passes it to Lake Erie, and so on. This creates a massive, coupled system of equations that is vital for predicting catastrophic flooding across the entire region.
___
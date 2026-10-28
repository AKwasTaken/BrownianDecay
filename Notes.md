## 7.6 Biological Waste Treatment, Hydrology and Indoor Air Quality

**Biological Waste Treatment** The continuous flow cultivation of microorganisms for the biological treatment of urban waste water has attracted widespread attention. Most of the mathematical models of microbial growth in current usage trace back to a model of Monod for single strain bacterial growth in a single substrate. Harris proposed a model of an anaerobic digester of sewage sludge based on a perfect mixing model of a continuously stirred tank reactor with noisy fluctuations in the concentrations and feedrates, and with flow equalization prior to the reaction being used to stabilize the process. He considered five state variables $(X^1, X^2, X^3, X^4, X^5)$ with $X^1$ the equalization tank volume, $X^2$ the incoming sludge flow rate, $X^3$ the incoming substrate (pollutant) concentration, $X^4$ the substrate concentration leaving the equalization tank for the reactor tank, and $X^5$ the substrate concentration leaving the reactor tank.

To account for random fluctuations in the influent flow rate and pollutant concentration due to rainfall and industrial and domestic waste discharge, Harris used the pair of Ito stochastic differential equations

$$
(6.1) \qquad dX_t^2 = \left(\alpha_1 + 1 - \beta_1 \frac{X_t^2}{X_t^1}\right) dt + \sqrt{2\gamma_1 \frac{X_t^2}{X_t^1}} dW_t^1
$$

$$
(6.2) \qquad dX_t^3 = \frac{(\alpha_2 + 1 - \beta_2 X_t^3)}{X_t^1} dt + \sqrt{2\gamma_2 \frac{X_t^3}{X_t^1}} dW_t^2
$$

to model the turbulent flow and concentration pattern. Here the $\alpha_i, \beta_i$ and $\gamma_i$ are positive parameters characterizing the noisy fluctuations. Then, using mass balance for the equalization tank, he obtained a pair of ordinary differential equations for the equalization tank volume and the substrate concentration leaving the equalization tank:

$$
(6.3) \qquad dX_t^1 = \left(\frac{\alpha_1 + 1}{\beta_1} - X_t^2\right) dt
$$

$$
(6.4) \qquad dX_t^4 = (X_t^3 - X_t^4) \frac{X_t^2}{X_t^1} dt.
$$

Finally, the mass balance for the reactor tank gave an ordinary differential equation for the substrate concentration leaving the reactor tank:

$$
(6.5) \qquad dX_t^5 = \left(a + b \left(X_t^4 - X_t^5\right)\right) dt
$$

where $a$ and $b$ are factors depending on the growth rate of the micro-organisms, the yield constant and the noise parameters.

The five equations (6.1)-(6.5) form a nonlinear 5-dimensional vector stochastic differential equation with no noise terms appearing explicitly in three of the five component equations, so the corresponding rows in the diffusion matrix have zero entries. Harris presented a numerical method for this system, and used it to plot the time evolution of the expected gain in concentration of the biological oxygen demanding states, that is $E(X_t^5/X_t^3)$ against time $t$ for various values of inputs and parameters.

**Hydrology** The outflow $x$ from a lake or reservoir depends on the precipitation $P$ over its catchment area, the evaporation from the water surface, the seepage flow and other factors which can be combined and represented abstractly by a parameter $\Phi$. Hydrologists have modelled the outflow from a single reservoir by a first order ordinary differential equation

$$
(6.6) \qquad \dot{x} = F(x)
$$

where the function on the right is given by

$$
(6.7) \qquad F(x) = AP x^\alpha - AP\Phi(1 - \psi(x)) x^\alpha \\ \qquad \qquad - AE\psi(x)x^\alpha - Ax^\beta - Bx^\gamma.
$$

The last term in (6.7) is the seepage flow and $\psi(x) = C + Dx^\delta$ is determined by the ratio of the surface areas of the catchment region and the reservoir. The constants $A$, $B$, $C$ and $D$ are nonnegative, whereas the indices $\alpha$, $\beta$, $\gamma$ and $\delta$ may also be negative.

Since the precipitation $P$, evaporation $E$ and catchment abstractions $\Phi$ are subject to fluctuations from many sources, it is reasonable to consider them to be random variables. Unny and Karmeshu argued that they could be approximated by independent Gaussian white noise fluctuations about positive mean values, namely

$$
P + \sigma_1 \xi_t^1, \quad E + \sigma_2 \xi_t^2 \quad \text{and} \quad \Phi + \sigma_3 \xi_t^3.
$$

They substituted these into (4.2), using $P\Phi + P\sigma_3\xi_t^3 + \Phi\sigma_1\xi_t^1$ for the product of $P + \sigma_1\xi_t^1$ and $\Phi + \sigma_3\xi_t^3$ in the second term, and obtained the Stratonovich stochastic differential equation

$$
(6.8) \qquad dX_t = F(X_t)dt + \sigma_1(1 - \Phi + A\Phi\psi(X_t)(X_t)^\alpha) \circ dW_t^1 \\ \qquad \qquad - \sigma_2 A \psi(X_t) \circ dW_t^2 \\ \qquad \qquad - \sigma_3 A P (1 - \psi(X_t))(X_t)^\alpha \circ dW_t^3
$$

as a stochastic counterpart of (6.6). Here $W_t^1$, $W_t^2$ and $W_t^3$ are independent standard Wiener processes. Unny and Karmeshu simulated numerically some sample paths of various solutions of (6.8). They also considered a cascade of lakes or reservoirs, such as the Great Lakes of North America, with the output of one being an additional input into the next. This resulted in a system of Stratonovich equations of the form (6.8) with additional coupling terms.

**Air Quality** Air conditioning units are used not only to regulate the temperature within a room, but also to filter the air to remove pollutants. Traditionally determinsitic models, based on mass and heat balances, have been used to help calibrate such units. Restricting attention here to the air quality in a single room, an established deterministic model for the concentration of a pollutant such as $CO_2$ is

$$
(6.9) \qquad v \dot{C}^i = -k \left(q^0 + f^1 q^1 + q^2\right) C^i + S - R \\ \qquad \qquad + k \left(q^0(1 - f^0) + q^2\right)C^0
$$

where $C^i$ is the indoor concentration and $C^0$ the outdoor concentration; $q$ is the volumetric flow rate for the make-up air $q^0$, the recirculation $q^1$, the infiltration $q^2$, the exfiltration $q^3$, and the exhaust $q^4$; $f$ is the filter efficiency for the makeup $f^0$ and the recirculation $f^1$ air; $v$ is the volume of the room; $S$ is the indoor source emission rate; $R$ is the indoor sink removal rate; and $k$ is a factor which accounts for the inefficiency of the mixing.

In reality, many of the above parameters cannot be predicted with certainty due to random fluctuations, which can be quite large. To obtain a more realistic model, Haghighat, Fazio and Unny replaced $S$, $C^0$, $q^0$ and $q^2$ by random fluctuations about deterministic values $\bar{S}$, $\bar{C}^0$, $\bar{q}^0$ and $\bar{q}^2$ (which may vary in time), that is by

$$
S_t = \bar{S} + \xi_t^{(1)}, \qquad C_t^0 = \bar{C}^0 + \xi_t^{(2)}, \\ q_t^0 = \bar{q}^0 + \xi_t^{(3)}, \qquad q_t^2 = \bar{q}^2 + \xi_t^{(4)},
$$

where $\xi_t^{(1)}, \xi_t^{(2)}, \xi_t^{(3)}$ and $\xi_t^{(4)}$ are noise terms. They substituted these into (6.9), neglecting the product $\xi_t^{(2)}\xi_t^{(3)}$, to obtain

$$
(6.10) \qquad v dC_t^i = g(C^i, t) dt + H(C^i, t) d\xi_t
$$

where $\xi_t = (\xi_t^{(1)}, \xi_t^{(2)}, \xi_t^{(3)}, \xi_t^{(4)})$, $g$ is the scalar function

$$
g(C^i, t) = -k \left(\bar{q}^0 + f^1 q^1 + \bar{q}^2\right) C^i + \bar{S} - R \\ + k \left(\bar{q}^0(1 - f^0) + \bar{q}^2\right)\bar{C}^0,
$$

and $H$ is the (row) vector function with

$$
H(C^i, t)^T = \begin{pmatrix} 1 \\ k \left(\bar{q}^0 (1 - f^0) + \bar{q}^2\right) \\ k(\bar{q}^0(1 - f^0) - C^i)\bar{C}^0 \\ k(\bar{C}^0 - C^i) \end{pmatrix}.
$$

Using independent Wiener processes for the noise terms, they interpreted (6.10) as an Ito stochastic differential equation. They also considered coupled systems of such equations to model the pollutant concentrations in several different rooms.

Haghighat, Chandrashekar and Unny proposed a model, with a similar mathematical structure, of the thermal variations within a building subject to external temperature fluctuations.
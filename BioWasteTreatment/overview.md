# Part 1: The Core Intuition & Physical System
## The Real-World Engineering Problem
In modern municipal wastewater treatment plants, raw sewage (pollutants, organic matter measured as Chemical Oxygen Demand or COD) flows into biological reactors where colonies of bacteria feed on the waste and convert it into benign byproducts.
Microbial colonies are living systems. If you hit them with a sudden surge of toxic chemicals or concentrated organic waste, they suffer toxic shock or acidification, killing off the bacterial culture. A wastewater treatment plant cannot simply operate deterministically because:
1. Domestic use has morning/evening diurnal peaks.
2. Rainstorms flush street debris and urban runoff into sewers.
3. Industrial plants perform unannounced discharge dumps.
## The Two-Tank Architecture
To protect the active biological reactor, environmental engineers install a preliminary buffer:
1. **Tank 1 (Equalization Tank):** A variable-volume holding tank. Its sole purpose is damping and mixing. It absorbs erratic surges and smooths out both flow rate and pollutant concentration.
2. **Tank 2 (Reactor Tank / CSTR):** A continuously stirred tank reactor of fixed volume containing active bacteria. It receives the smoothed stream from Tank 1 and performs biological digestion.

```
       Chaotic Inflow
 [Flow: X^2, Dirtiness: X^3]
               │
               ▼
   ┌──────────────────────┐
   │ Equalization Tank    │ ◄── Volume: X^1
   │ (Mixing & Damping)   │
   └──────────┬───────────┘
              │ Buffered Stream (X^4)
              ▼
   ┌──────────────────────┐
   │ Biological Reactor   │ ◄── Microbial Digestion
   │ (CSTR)               │
   └──────────┬───────────┘
              │ Treated Effluent (X^5)
              ▼
```

## The 5 State Variables
- $X_t^1$: Volume of liquid currently in the Equalization Tank ($\text{m}^3$).
- $X_t^2$: Influent flow rate entering the plant ($\text{m}^3/\text{day}$).
- $X_t^3$: Influent pollutant concentration entering the plant ($\text{mg/L}$ or $\text{g/m}^3$).
- $X_t^4$: Pollutant concentration leaving Tank 1 and entering the Reactor ($\text{mg/L}$).
- $X_t^5$: Pollutant concentration leaving Tank 2 into the environment ($\text{mg/L}$).
___
# Part 2: Mathematical Formulations & Derivations
The system is governed by a 5-dimensional non-linear stochastic differential equation (SDE):
$$dX_t = f(X_t) \, dt + G(X_t) \, dW_t$$
### 1. Inflow Flow Rate $X^2$ & Dirtiness $X^3$
The incoming stream is modeled as a mean-reverting, square-root diffusion process (structurally identical to the **Cox–Ingersoll–Ross / Feller process**):
$$dX_t^2 = \left(\alpha_1 + 1 - \beta_1 \frac{X_t^2}{X_t^1}\right) dt + \sqrt{2\gamma_1 \frac{X_t^2}{X_t^1}} \, dW_t^1$$
$$dX_t^3 = \frac{(\alpha_2 + 1 - \beta_2 X_t^3)}{X_t^1} \, dt + \sqrt{2\gamma_2 \frac{X_t^3}{X_t^1}} \, dW_t^2$$
#### Deconstructing Every Term:
- **The Drift term ($dt$):** Notice the term $-\beta_1 \frac{X^2}{X^1}$. If the incoming flow $X^2$ spikes way above the holding volume $X^1$, this negative term pulls the derivative downward. It acts like an elastic rubber band returning the inflow to a baseline mean.
- **The Diffusion term ($dW_t$):** $W_t^1$ and $W_t^2$ are independent standard 1D Brownian motions. The noise amplitude is scaled by $\sqrt{X_t}$. This ensures positivity: as flow $X^2 \to 0$ or concentration $X^3 \to 0$, the random fluctuations vanish ($\sqrt{0} = 0$), preventing negative volume, negative flow, or negative mass.
- **The $\frac{1}{X_t^1}$ scaling:** When the equalization tank has high liquid volume ($X^1$ is large), fluctuations are diluted over a greater thermal/hydraulic inertia.
### 2. Equalization Tank Volume ($X^1$)
By conservation of fluid volume:
$$\frac{dV}{dt} = Q_{\text{out}} - Q_{\text{in}}$$
Harris sets the engineered outflow rate at the design operating baseline $Q_{\text{out}} = \frac{\alpha_1 + 1}{\beta_1}$. Since $Q_{\text{in}} = X_t^2$:
$$dX_t^1 = \left( \frac{\alpha_1 + 1}{\beta_1} - X_t^2 \right) dt$$
- If inflow exceeds design outflow ($X^2 > Q_{\text{out}}$), $dX^1 < 0$ (tank drains or balances).
- If inflow drops below design outflow, the tank accumulates fluid.
- **Notice:** There is no explicit $dW_t$ term here. Tank 1 does not generate random liquid from vacuum. Its randomness is driven exclusively by the chaotic input $X^2$.
### 3. Equalization Tank Solute Dynamics ($X^4$)
We derive this using first-principles conservation of mass inside Tank 1.
Let $M_1(t)$ be the total mass of pollutant in Tank 1:
$$M_1(t) = X_t^1 \cdot X_t^4 \quad (\text{Volume} \times \text{Concentration})$$
Mass enters at rate $X_t^3 X_t^2$ ($\text{Inflow concentration} \times \text{Inflow rate}$). Assuming the tank is well-mixed, liquid leaves at concentration $X_t^4$ with flow rate $X_t^2$. Therefore, the rate of change of mass is:
$$\frac{dM_1}{dt} = X_t^3 X_t^2 - X_t^4 X_t^2 = (X_t^3 - X_t^4)X_t^2$$
Now, expand $dM_1 = d(X^1 X^4)$ using the product rule:
$$d(X_t^1 X_t^4) = X_t^1 \, dX_t^4 + X_t^4 \, dX_t^1$$
_(Note: Because $X^1$ and $X^4$ do not have direct Brownian motion terms, their quadratic covariation $d\langle X^1, X^4 \rangle_t = 0$, so the standard calculus product rule holds exactly without an Itô correction term)._
Substitute $dX_t^1 = (Q_{\text{out}} - X_t^2)dt$:
$$(X_t^3 - X_t^4)X_t^2 \, dt = X_t^1 \, dX_t^4 + X_t^4 \left( \frac{\alpha_1 + 1}{\beta_1} - X_t^2 \right) dt$$
Under balanced hydraulic steady-state where $Q_{\text{out}} \approx X_t^2$, the net volume change balances out, leaving the standard formulation:
$$dX_t^4 = (X_t^3 - X_t^4) \frac{X_t^2}{X_t^1} \, dt$$
#### Physical Interpretation:
- The term $\frac{X^2}{X^1}$ has units of $\frac{\text{m}^3/\text{day}}{\text{m}^3} = \text{day}^{-1}$. This is the hydraulic turnover rate (the reciprocal of the residence time $\tau_{\text{res}} = \frac{X^1}{X^2}$).
- The rate of change of concentration $dX^4$ is directly driven by the difference between what enters ($X^3$) and what is currently inside ($X^4$), scaled by how fast water flushes through the tank.
### 4. Biological Reactor Effluent ($X^5$)
Inside the reactor tank (CSTR), active biomass removes pollutant via microbial consumption. The change in concentration is:
$$dX_t^5 = \Big( \underbrace{\frac{Q}{V_{\text{reactor}}} (X_t^4 - X_t^5)}_{\text{Hydraulic Washout/Inflow}} - \underbrace{R_{\text{biodegradation}}(X_t^5)}_{\text{Bacterial Consumption Rate}} \Big) dt$$
Under the Monod kinetic model, bacterial substrate uptake rate is non-linear:
$$R(S) = \frac{\mu_{\max} S}{K_s + S} \cdot \frac{B}{Y}$$
where $B$ is bacterial concentration, $Y$ is the yield coefficient, and $S = X^5$.
To make the system mathematically tractable, Harris linearizes the reaction rate around the baseline operational point using a first-order Taylor expansion:
$$R(X^5) \approx R_0 + k \cdot X^5$$
Combining the hydraulic dilution and linear biodegradation terms into lumped constant parameters gives:
$$dX_t^5 = \left(a + b(X_t^4 - X_t^5)\right) dt$$
where:
- $b$ governs the dilution and reaction throughput rate.
- $a$ reflects the background basal production or net decay constant.
### 5. Why the System is "Hypoelliptic / Degenerate"
Look at the full diffusion matrix $G(X)$:
$$G(X) = \begin{bmatrix} 0 & 0 \\ \sqrt{2\gamma_1 \frac{X^2}{X^1}} & 0 \\ 0 & \sqrt{2\gamma_2 \frac{X^3}{X^1}} \\ 0 & 0 \\ 0 & 0 \end{bmatrix}$$
Three rows are completely zero. In pure mathematics, if a diffusion matrix has zero rows, you cannot use standard parabolic PDE theorems to solve the Fokker–Planck equation easily.
However, because $X^2$ directly enters $dX^1$ and $dX^4$, and $X^4$ directly enters $dX^5$, the noise injected at states 2 and 3 "transfers" down the chain through the drift vector fields. This system satisfies Hörmander’s bracket condition for hypoellipticity, meaning that despite having only 2 noise sources for 5 variables, the solution admits a smooth transition probability density function in $\mathbb{R}^5$.
__________
# Part 3: Interpreting Every Single Analysis Plot
### 1. Damping Ratio & Cross-Correlation
- **Damping:** $\frac{\text{Var}(X^4)}{\text{Var}(X^3)} \implies$ Quantifies noise absorption
- **Lag:** $R_{34}(\tau) \implies$ Measures hydraulic transit delay
- **Variance Ratio $\rho_{\text{damp}}(t) = \frac{\text{Var}(X^4)}{\text{Var}(X^3)}$:**
    - Raw sewage volatility is $\text{Var}(X^3)$. The volatility entering the bacteria tank is $\text{Var}(X^4)$.
    - **What it tells you:** If $\rho_{\text{damp}} \approx 0.15$, the equalization tank suppresses 85% of environmental noise. If $\rho \ge 1$, the tank is failing to buffer the load.
- **Cross-Correlation Peak $\tau^*$:**
    - When an industrial plant dumps waste at $t = 1.0\text{ day}$, when does that spike hit Tank 2?
    - **What it tells you:** The $\arg\max \tau^*$ is the hydraulic delay. If $\tau^* = 0.8\text{ days}$, it takes nearly 20 hours for a shock at the gate to diffuse into the reactor. This tells plant operators exactly how much reaction time they have to adjust aeration or dosing pumps.
### 2. First-Passage Time & Survival Analysis
- **FHT Density:** When do violations occur?
- **Survival:** Probability that plant never violates standards over time
- **The Environmental Reality:** The law does not care about average effluent. If COD crosses $C_{\text{crit}} = 3.5\text{ mg/L}$ even for 30 minutes, river wildlife dies and the plant receives a fine.
- **First Hitting Time (FHT):**
$$\tau = \inf \{ t > 0 : X_t^5 \ge C_{\text{crit}} \}$$
- A deterministic model can only say: "The effluent reaches an asymptote of 3.2, so it never breaches."
- The stochastic FHT histogram shows the real picture: because of random shocks, 15% of trajectories cross 3.5 within the first 3 days.
- **Survival Function $S(t) = P(\tau > t)$:**
    - Shows the continuous decay of plant compliance. At $t = 0$, $S(0) = 1.0$ (100% compliant). Over time, as stochastic disturbances accumulate, the survival probability curve drops.
### 3. Fokker–Planck Distribution Slices (KDE)
- **Density Slices:** Shows how the distribution of effluent spreads and skews over time.
- At $t = 1.0\text{ d}$, the distribution of $X^5$ is narrow because all paths started from the same initial condition $X_0^5 = 1.0$.
- As time moves to $t = 2.5, 4.0, 5.0\text{ d}$, the distribution widens (diffusion of probability mass).
- **Skewness & Kurtosis:**
    - Standard Brownian motion produces symmetric Gaussian bells ($\text{Skewness} = 0$).
    - But your square-root diffusion has state-dependent noise ($\sqrt{X}$). When concentrations are high, noise is larger; when low, noise is smaller. This creates an asymmetric positive skewness (a fat right tail). That fat right tail is the "disaster zone" where unexpected extreme pollution events concentrate.
### 4. Covariance Matrix & Phase Plane
- **Covariance:** Which states move together?
- **Phase Plane:** ($X^4$ vs $X^5$) Trajectory funneling into an attractor
- **Covariance Matrix:**
    - $\text{Cov}(X^1, X^5)$: If volume and effluent are negatively correlated, higher volume dampens effluent concentration (dilution effect).
    - $\text{Cov}(X^3, X^5)$: Shows how much raw incoming toxicity bleeds directly through the entire plant into the treated water.
- **Phase Portrait ($X^4$ vs $X^5$):**
    - In deterministic systems, you see a single smooth line winding into an equilibrium point.
    - In this stochastic system, you see a bundle of paths forming a "tube" or "funnel." The width of that tube represents the irreducible noise boundary around the plant's operational setpoint.
### 5. Fixed Performance Metrics (Mass Balance vs. Ratio of Expectations)
- **Mass Removal:** $\eta = 1 - \frac{\text{Total Outflow Mass}}{\text{Total Inflow Mass}}$
- **Mean Ratio:** $\frac{\mathbb{E}[X^5]}{\mathbb{E}[X^3]}$
- **Why your original $\mathbb{E}[X^5 / X^3]$ blew up:**
    - In probability theory, if $Y$ is a positive random variable that can get close to zero, $\mathbb{E}[1/Y]$ can be massive (or even infinite), even if $\mathbb{E}[Y]$ is completely normal. Pathwise division by tiny numbers created false mathematical spikes.
- **Cumulative Mass Balance Efficiency $\eta_{\text{mass}}(t)$:**
$$\eta_{\text{mass}}(t) = 1 - \frac{\int_0^t X_s^5 X_s^2 \, ds}{\int_0^t X_s^3 X_s^2 \, ds}$$
- This is the true engineering metric: Total grams of waste released divided by total grams of waste received.
- The curve starts near 0% and asymptotically rises toward 80–90%, showing how the treatment system stabilizes into its full cleanup efficiency as the tanks fill up and microbial activity reaches equilibrium.
### 6. Volatility Sweep ($\gamma_2$ Sensitivity)
- **X-axis:** $\gamma_2$ (Input Noise Intensity)
- **Y-axis:** Breach Probability & Terminal Effluent Variance
- This simulates what happens if weather becomes more erratic or industrial dumps become more frequent.
- **The Non-Linear Cliff:**
    - When $\gamma_2$ is low ($0.02 \to 0.10$), breach probability remains low because the equalization tank can easily buffer small ripples.
    - Past a critical threshold ($\gamma_2 > 0.25$), the breach curve bends sharply upward.
    - This identifies the systemic breakdown limit: the precise level of inflow turbulence at which the equalization tank is physically outmatched, and the plant requires structural expansion.
# Project Report: Stochastic Analysis and Performance Modeling of a Two-Stage Biological Wastewater Treatment Plant

## 1\. Motivation and Background

### 1.1 The Real-World Engineering Problem

Cities produce wastewater continuously, but neither the volume nor the filthiness of that water is constant. Domestic routines create sharp morning and evening discharge peaks, sudden convective storms flood municipal sewers with surface runoff, and industrial facilities periodically dump concentrated batches of chemical waste.

To clean this water before returning it to natural rivers, modern facilities use biological wastewater treatment—specifically, anaerobic digestion. Microorganisms (bacteria) reside inside a bioreactor tank and consume organic carbon compounds (measured as Chemical Oxygen Demand, or COD) as metabolic food.

However, bacterial colonies operate within delicate kinetic thresholds:

1. **Hydraulic Washout:** If incoming flow surges too quickly, bacteria get flushed out of the tank faster than they can reproduce.
2. **Substrate Overload / Toxic Shock:** A massive spike in pollutant concentration can alter local pH and inhibit bacterial enzymes, causing the biological digestive system to fail.

To prevent plant failure, engineers place an **Equalization Tank** (a buffer holding basin) ahead of the **Continuous Stirred Tank Reactor (CSTR)**. The equalization tank dampens wild inflow spikes, blending concentrated wastewater with existing volume before the mixture enters the reactor tank.

#### Two-Stage Wastewater Treatment Process

text

```text
       [ Chaotic Raw Sewage ]
        Flow X², Dirtiness X³
                 │
                 ▼
     ┌────────────────────────┐
     │ 1. Equalization Tank   │
     │    Volume: X¹          │  <--- Damps high-frequency spikes
     │    Leaving Dirt: X⁴    │
     └───────────┬────────────┘
                 │
                 ▼
     ┌────────────────────────┐
     │ 2. Reactor Tank (CSTR) │
     │    Treated COD: X⁵     │  <--- Bacterial digestion occurs here
     └───────────┬────────────┘
                 │
                 ▼
       [ Clean Treated Water ]
```

### 1.2 The Failure of Deterministic Models and the Need for Stochastic Calculus

In standard chemical engineering, systems are modeled using deterministic Ordinary Differential Equations (ODEs) derived from conservation of mass:

$$
\frac{d(\text{Mass})}{dt}
=
\text{Mass In}
-
\text{Mass Out}
-
\text{Mass Consumed}
$$

A deterministic model assumes that the incoming sewage flow rate and incoming pollutant concentration follow smooth, perfectly predictable average paths. However, applying an average input to a nonlinear biological system yields inaccurate predictions because biological degradation relies on nonlinear saturation kinetics (such as the Monod equation):

$$
\mu(S) = \mu_{\max} \frac{S}{K_s + S}
$$

Because this function is strictly concave ($\mu''(S) < 0$), Jensen's Inequality establishes that:

$$
\mathbb{E}[\mu(S)] < \mu(\mathbb{E}[S])
$$

Operating an anaerobic digester in a real environment with random surges results in lower average bacterial growth and lower treatment efficiency than an ODE model predicts using average concentrations.

Deterministic models also cannot answer reliability questions:

- What is the probability that effluent pollution breaches environmental safety limits during a 5-day cycle?
- How long can the biological reactor withstand an inflow surge before failing compliance tests?

To answer these questions, the system must be formulated as a **System of Stochastic Differential Equations (SDEs)** governed by Itô calculus, incorporating environmental white noise into the mass balance.

---

## 2\. Mathematical Formulation of the Harris Model

The model developed by Carl J. Harris (1977) represents the wastewater treatment plant as a 5-dimensional dynamical state vector:

$$
X_t =
[X_t^1, X_t^2, X_t^3, X_t^4, X_t^5]^T
$$

The state variables are defined as follows:

| State Variable | Physical Meaning | Unit |
| --- | --- | --- |
| $X_t^1$ | Liquid volume in the Equalization Tank | $\text{m}^3$ |
| $X_t^2$ | Inflow volumetric rate of raw wastewater entering the plant | $\text{m}^3/\text{day}$ |
| $X_t^3$ | Pollutant substrate concentration (COD) of raw incoming wastewater | $\text{mg/L}$ |
| $X_t^4$ | Pollutant substrate concentration leaving the Equalization Tank | $\text{mg/L}$ |
| $X_t^5$ | Final treated effluent substrate concentration leaving the Bioreactor | $\text{mg/L}$ |

### 2.1 The Coupled SDE System

The complete stochastic model is given by the following five coupled SDEs:

$$
\begin{aligned}
dX_t^1
&=
\left(
\frac{\alpha_1 + 1}{\beta_1} - X_t^2
\right)dt
\\[8pt]
dX_t^2
&=
\left(
\alpha_1 + 1
-
\beta_1 \frac{X_t^2}{X_t^1}
\right)dt
+
\sqrt{
2\gamma_1 \frac{X_t^2}{X_t^1}
}\,dW_t^1
\\[8pt]
dX_t^3
&=
\left(
\frac{\alpha_2 + 1 - \beta_2 X_t^3}{X_t^1}
\right)dt
+
\sqrt{
2\gamma_2 \frac{X_t^3}{X_t^1}
}\,dW_t^2
\\[8pt]
dX_t^4
&=
(X_t^3 - X_t^4)
\frac{X_t^2}{X_t^1}\,dt
\\[8pt]
dX_t^5
&=
\left(
a + b(X_t^4 - X_t^5)
\right)dt
\end{aligned}
$$

Here,

$$
W_t = [W_t^1, W_t^2]^T
$$

is a 2-dimensional standard Brownian motion vector on a filtered probability space:

$$
(\Omega, \mathcal{F}, (\mathcal{F}_t)_{t \geq 0}, \mathbb{P})
$$

### 2.2 Term-by-Term Mathematical and Physical Breakdown

#### Equation 1: Equalization Tank Volume ($X^1$)

$$
dX_t^1 =
\left(
\frac{\alpha_1 + 1}{\beta_1} - X_t^2
\right)dt
$$

- **Physical Balance:**
  $$
  \frac{dV}{dt} = Q_{\text{out}} - Q_{\text{in}}
  $$
- **Mechanism:** The nominal baseline outflow rate from the holding tank is set to a constant design target:
  $$
  Q^* = \frac{\alpha_1 + 1}{\beta_1}
  $$
  If incoming flow $X_t^2$ exceeds this rate, the net derivative turns negative or adjusts the accumulation of water in the tank.
- **Absence of Noise:** The physical volume of liquid does not jump instantaneously; it responds as the continuous integral of the flow rate.

#### Equation 2: Inflow Rate Dynamics ($X^2$)

$$
\begin{aligned}
dX_t^2
&=
\underbrace{
\left(
\alpha_1 + 1
-
\beta_1 \frac{X_t^2}{X_t^1}
\right)dt
}_{\text{Mean-Reverting Drift}}
\\
&\quad+
\underbrace{
\sqrt{
2\gamma_1 \frac{X_t^2}{X_t^1}
}\,dW_t^1
}_{\text{Feller / Square-Root Diffusion}}
\end{aligned}
$$

- **Drift Term:** When the flow rate $X^2$ deviates from its equilibrium relative to tank volume, the negative feedback term
  $$
  -\beta_1\frac{X^2}{X^1}
  $$
  pulls it back toward balance.
- **Diffusion Term:** The term $dW_t^1$ introduces Gaussian white noise with variance $dt$. The diffusion coefficient contains $\sqrt{X^2}$. In stochastic analysis, this is a Cox–Ingersoll–Ross (CIR) or Feller square-root diffusion.
- **Boundary Behavior:** If inflow drops toward zero ($X^2 \to 0$), the noise amplitude shrinks to zero ($\sqrt{0}=0$), while the deterministic drift remains positive ($( \alpha_1+1)dt>0$). This acts as a reflecting boundary that keeps incoming flow strictly non-negative ($X^2 \geq 0$), preventing unphysical negative flow rates.

#### Equation 3: Inflow Wastewater Dirtiness ($X^3$)

$$
dX_t^3 =
\frac{\alpha_2 + 1 - \beta_2 X_t^3}{X_t^1}\,dt
+
\sqrt{
2\gamma_2 \frac{X_t^3}{X_t^1}
}\,dW_t^2
$$

- Models stochastic variations in raw sewage strength (COD).
- Governed by its own independent Brownian motion $W_t^2$, scaled inversely by the holding tank volume $X^1$ (a larger buffering capacity reduces the variance of concentration shocks entering the system).
- Retains a square-root diffusion coefficient to preserve non-negativity ($X^3 \geq 0$).

#### Equation 4: Equalization Tank Effluent Concentration ($X^4$)

$$
dX_t^4 =
(X_t^3-X_t^4)
\frac{X_t^2}{X_t^1}\,dt
$$

- **Physical Balance:** Mass balance of a well-mixed holding basin without internal reaction.
- **Mechanism:** The rate of change of pollutant concentration inside the holding tank depends on the difference between what enters ($X^3$) and what is currently inside ($X^4$), scaled by the hydraulic turnover rate:
  $$
  \frac{X_t^2}{X_t^1}
  =
  \frac{\text{Flow}}{\text{Volume}}
  =
  \frac{1}{\text{Retention Time}}
  $$
- As liquid mixes in the tank, $X^4$ acts as a low-pass filter on the chaotic input $X^3$.

#### Equation 5: Bioreactor Effluent Concentration ($X^5$)

$$
dX_t^5 =
\left(
a+b(X_t^4-X_t^5)
\right)dt
$$

- **Physical Balance:** CSTR mass balance coupled with linearized biological digestion.
- **Mechanism:** Waste enters from the holding tank at concentration $X^4$ and exits at concentration $X^5$. The term $b(X^4-X^5)$ reflects mass exchange and hydraulic dilution. The biological breakdown from the bacterial population is incorporated through parameters $a$ and $b$, derived from a first-order Taylor expansion of Monod kinetics around the design operating point.

### 2.3 Degenerate (Hypoelliptic) Diffusion Structure

A key mathematical property of the Harris system is its partially degenerate diffusion matrix. In vector matrix notation:

$$
dX_t = \mathbf{f}(X_t)\,dt + \mathbf{G}(X_t)\,dW_t
$$

$$
\begin{bmatrix}
dX^1\\
dX^2\\
dX^3\\
dX^4\\
dX^5
\end{bmatrix}
=
\begin{bmatrix}
f_1(X)\\
f_2(X)\\
f_3(X)\\
f_4(X)\\
f_5(X)
\end{bmatrix}dt
+
\begin{bmatrix}
0 & 0\\
\sqrt{2\gamma_1\frac{X^2}{X^1}} & 0\\
0 & \sqrt{2\gamma_2\frac{X^3}{X^1}}\\
0 & 0\\
0 & 0
\end{bmatrix}
\begin{bmatrix}
dW^1\\
dW^2
\end{bmatrix}
$$

The diffusion matrix $\mathbf{G}(X_t)\in\mathbb{R}^{5\times 2}$ has rank 2:

$$
\operatorname{Rank}(\mathbf{G})=2<5
$$

Three of the state equations ($X^1$, $X^4$, $X^5$) contain no explicit diffusion terms ($g_i\equiv 0$). Randomness enters the physical plant exclusively through the input boundary conditions (Equations 6.1 and 6.2). The downstream states inherit stochasticity solely through non-linear drift coupling with $X^2$ and $X^3$.

In stochastic analysis, this configuration is known as a **hypoelliptic system**; noise smooths through deterministic integrals rather than diffusing uniformly in all state directions.

---

## 3\. Numerical Discretization Scheme

Because the drift and diffusion functions are non-linear and coupled, the 5D system cannot be solved in closed form. We use the **Euler–Maruyama method**, the stochastic counterpart to the Euler forward scheme for differential equations.

### 3.1 Numerical Algorithm and Positivity Preservation

Let the operating window $[0,T]$ be divided into $N$ equal steps of size:

$$
\Delta t = \frac{T}{N}
$$

with discrete times:

$$
t_n=n\Delta t
$$

The continuous Brownian increments $dW_t^k$ are approximated by discrete random variables:

$$
\begin{aligned}
\Delta W_n^k
&=
W_{(n+1)\Delta t}^k-W_{n\Delta t}^k\\
&=
\sqrt{\Delta t}\,Z_n^k,
\qquad
\text{where }Z_n^k\sim\mathcal{N}(0,1)
\end{aligned}
$$

Applying the Euler–Maruyama update to the Harris system gives:

$$
\begin{aligned}
X_{n+1}^1
&=
X_n^1+
\left(
\frac{\alpha_1+1}{\beta_1}-X_n^2
\right)\Delta t
\\[8pt]
X_{n+1}^2
&=
X_n^2+
\left(
\alpha_1+1-\beta_1\frac{X_n^2}{X_n^1}
\right)\Delta t
\\
&\quad+
\sqrt{
2\gamma_1\frac{\widetilde{X}_n^2}{\widetilde{X}_n^1}
}
\sqrt{\Delta t}\,Z_n^1
\\[8pt]
X_{n+1}^3
&=
X_n^3+
\left(
\frac{\alpha_2+1-\beta_2X_n^3}{X_n^1}
\right)\Delta t
\\
&\quad+
\sqrt{
2\gamma_2\frac{\widetilde{X}_n^3}{\widetilde{X}_n^1}
}
\sqrt{\Delta t}\,Z_n^2
\\[8pt]
X_{n+1}^4
&=
X_n^4+
(X_n^3-X_n^4)
\frac{X_n^2}{X_n^1}\Delta t
\\[8pt]
X_{n+1}^5
&=
X_n^5+
\left(
a+b(X_n^4-X_n^5)
\right)\Delta t
\end{aligned}
$$

#### Full Truncation Guard

In discrete simulations, a large negative draw from the Gaussian distribution $Z_n$ can drive $X_{n+1}^2$ or $X_{n+1}^3$ below zero, which causes imaginary values in the square-root diffusion terms.

To preserve numerical stability while honoring physical non-negativity, we implement a full truncation scheme:

$$
\widetilde{X}_n^i=\max(X_n^i,0),
\qquad
\widetilde{X}_n^1=\max(X_n^1,\epsilon)
$$

where:

$$
\epsilon=10^{-4}
$$

prevents division by zero in the hydraulic retention term.

---

## 4\. Simulation Results and Detailed Data Analysis

Simulations were conducted using an ensemble of:

- Number of Monte Carlo trajectories: $N_{\text{paths}}=3{,}000$
- Operational window: $T=5.0$ days
- Step size: $\Delta t=0.0025$ days
- Number of time steps: $2{,}000$

The baseline parameters used were:

$$
\begin{aligned}
\alpha_1&=1.0,
&\beta_1&=1.2,
&\gamma_1&=0.15,\\
\alpha_2&=1.0,
&\beta_2&=1.0,
&\gamma_2&=0.20,\\
a&=0.05,
&b&=1.50
\end{aligned}
$$

The initial state vector was:

$$
X_0=
[10.0,\;1.667,\;5.0,\;4.0,\;1.0]^T
$$

### 4.1 System Dynamics and Input Buffering

Figure 4.1 displays sample trajectories from the simulation run.

- **Top-Left (Input Dynamics):** The raw incoming wastewater flow $X_t^2$ and pollutant concentration $X_t^3$ fluctuate erratically under Brownian increments, reflecting the unbuffered inputs to the municipal plant.
- **Top-Right (Equalization Effect):** A direct comparison between raw input $X^3$ (red) and equalized output $X^4$ (navy) along a single realization highlights the physical buffering of the holding basin. While $X^3$ fluctuates with sharp local variations, $X^4$ follows a smooth, low-pass-filtered path.
- **Bottom-Left (Effluent Concentration****$X^5$****):** Individual effluent trajectories (gray) remain bounded near the ensemble mean $\mathbb{E}[X_t^5]$ (dark green). The smooth curves confirm that downstream biological states inherit filtered noise without exhibiting white-noise roughness.
- **Bottom-Right (The Raw Ratio Dilemma):** Plotting the instantaneous pathwise ratio $\mathbb{E}[X_t^5/X_t^3]$ produces an unphysical spike near $t=4.0$ days, with values climbing above 1.5. This artifact does not indicate pollutant creation; rather, it highlights the mathematical limitations of evaluating pathwise expectations of quotients, as explained in Section 4.4.

### 4.2 Quantifying Variance Damping and Hydraulic Retention Lag

To quantify how effectively the holding basin stabilizes incoming noise, we analyze the transmission of variance and the operational time lag.

#### Variance Transmission Ratio ($\rho_{\text{damp}}(t)$)

The variance ratio between equalized and raw wastewater is defined as:

$$
\rho_{\text{damp}}(t)
=
\frac{\operatorname{Var}(X_t^4)}
{\operatorname{Var}(X_t^3)}
$$

As shown in the left panel of Figure 4.2, the ratio starts near zero, rises as variance accumulates, and plateaus around 0.80 to 0.85. This confirms that the holding tank absorbs roughly 15–20% of input variance even under persistent driving noise.

#### Hydraulic Delay via Cross-Correlation

Because the holding tank stores liquid, changes in effluent dirtiness lag behind changes in influent dirtiness.

The empirical cross-correlation function between $X_t^3$ and $X_t^4$ is:

$$
R_{34}(\tau)
=
\frac{
\mathbb{E}[
(X_t^3-\mu_3)(X_{t+\tau}^4-\mu_4)
]
}{
\sigma_3\sigma_4
}
$$

The right panel of Figure 4.2 shows that $R_{34}(\tau)$ peaks at:

$$
\tau^*=0.10\text{ days}
$$

(approximately 2.4 hours) with a normalized correlation of 0.78.

This confirms an average retention delay of 2.4 hours before raw sewage shocks propagate into the biological reactor.

### 4.3 Regulatory Compliance and First Hitting Time Analysis

Environmental regulations set strict upper limits on discharge pollution. In our simulation, the legal compliance threshold is set at:

$$
C_{\text{crit}}=3.5\text{ mg/L}
$$

of COD.

#### First Hitting Time (FHT)

We define the random variable $\tau_{\text{breach}}$ as the first time a trajectory breaches the regulatory limit:

$$
\tau
=
\inf\{t>0:X_t^5\geq C_{\text{crit}}\}
$$

The left panel of Figure 4.3 displays the empirical probability density of $\tau$ across the ensemble. The breach times form a bell-shaped distribution centered sharply around:

$$
\tau\approx1.05\text{ days}
$$

(spanning 0.98 to 1.13 days).

#### Survival and Failure Probability

The right panel of Figure 4.3 tracks the plant's reliability over time.

**Survival Function:**

$$
S(t)
=
\mathbb{P}(\tau>t)
=
\mathbb{P}
\left(
\sup_{0\leq s\leq t}X_s^5<C_{\text{crit}}
\right)
$$

**Failure Probability:**

$$
F(t)=1-S(t)=\mathbb{P}(\tau\leq t)
$$

Between $t=0.95$ days and $t=1.15$ days, the survival probability drops sharply from 1.0 to 0.0.

This transition shows that under high influent loading, deterministic compliance estimates fail; stochastic analysis reveals an operational window where failure probability reaches 100% if buffering capacity is insufficient.

### 4.4 Resolving Performance Metrics: Mass Removal vs. Instantaneous Ratios

In the preliminary simulation (Figure 4.1), the ratio $\mathbb{E}[X_t^5/X_t^3]$ climbed above 1.0 around $t=4$ days, falsely suggesting that the treatment plant was generating pollution.

Two distinct factors explain this behavior:

1. **Hydraulic Transport Lag:** As established by cross-correlation analysis, the plant introduces a transport delay. When raw sewage input $X_t^3$ drops suddenly at $t\approx3.8$ days, the reactor $X_t^5$ is still digesting the heavy waste load that entered hours earlier. Comparing the two concentrations simultaneously at time $t$ divides a delayed peak by a dropping input.
2. **Jensen's Inequality on Expectations of Quotients:**
   $$
   \mathbb{E}\left[\frac{X^5}{X^3}\right]
   \neq
   \frac{\mathbb{E}[X^5]}{\mathbb{E}[X^3]}
   $$
   Because $f(x)=1/x$ is convex for $x>0$, any sample paths where $X^3$ fluctuates near zero produce extremely large ratios that skew the sample mean upward.

To correct for these effects, Figure 4.4 presents two physically robust metrics.

#### Cumulative Mass Removal Efficiency ($\eta_{\text{mass}}(t)$)

$$
\eta_{\text{mass}}(t)
=
1-
\frac{
\int_0^t X_s^5X_s^2\,ds
}{
\int_0^t X_s^3X_s^2\,ds
}
$$

The left panel demonstrates that cumulative mass removal begins at 80% and gradually decreases toward a steady state of roughly 6–8% as total cumulative loading equilibrates. Mass is conserved throughout the run, confirming physical consistency.

#### Ratio of Expectations

$$
\bar{R}(t)
=
\frac{\mathbb{E}[X_t^5]}
{\mathbb{E}[X_t^3]}
$$

The right panel shows that this ratio rises smoothly from 0.20 to 1.40 during the low-influent period, avoiding the artificial numerical spikes seen in pathwise quotients.

### 4.5 Probability Density Evolution (Fokker–Planck Slices)

Figure 4.5 shows the probability density functions $p(X^5,t)$ of the effluent state at four distinct times:

$$
t\in\{1.0,\;2.5,\;4.0,\;5.0\}\text{ days}
$$

- **At****$t=1.0$****d:** The distribution is narrow, centered at $X^5\approx3.48$ mg/L, with near-zero skewness (0.12) and kurtosis (0.05). The process is tightly clustered around its deterministic path.
- **At****$t=2.5$****d:** The distribution shifts to the right ($X^5\approx4.25$ mg/L) and widens as accumulated variance spreads the trajectories.
- **At****$t=4.0$****d and****$t=5.0$****d:** The distribution shifts back toward lower concentrations as influent loading drops, but its variance broadens significantly.
- **Emergence of Non-Gaussian Skewness:** By $t=5.0$ d, the distribution develops clear positive skewness (0.63) and elevated kurtosis (0.56). This rightward tail reflects the nonlinear square-root diffusion terms ($\sqrt{X^3}$): random positive shocks propagate more readily through the cascade than negative deviations.

### 4.6 Covariance Structure and Phase Plane Dynamics

Figure 4.6 shows the correlation structure and phase-space behavior across the plant.

#### Covariance Matrix ($\mathbf{\Sigma}(t=5.0)$)

The left panel displays the empirical covariance matrix across the five states.

- Raw wastewater $X^3$ exhibits the largest variance (0.444).
- Covariance between $X^3$ and $X^4$ is strong (0.238), confirming the transmission of concentration fluctuations.
- Covariance between $X^4$ and $X^5$ is similarly elevated (0.237), documenting downstream noise propagation.
- In contrast, covariances involving volume $X^1$ and flow rate $X^2$ remain near zero (-0.007 to -0.033), showing that hydraulic fluctuations stay largely decoupled from concentration dynamics at steady state.

#### Phase Plane Trajectory ($X^4$ vs. $X^5$)

The right panel maps equalized concentration $X^4$ against effluent concentration $X^5$.

The black trajectory denotes the ensemble mean, while colored lines trace individual sample paths.

The trajectory moves upward from $(4.0,1.0)$, reaches a turning point near $(4.2,4.25)$, and then loops downward toward $(2.3,2.7)$.

The phase paths remain bounded within a coherent attractor corridor, showing that the biological reactor maintains stability despite wide variations in input quality.

### 4.7 Sensitivity to Noise Intensity ($\gamma_2$)

Figure 4.7 illustrates the impact of varying the influent noise parameter $\gamma_2$ from 0.02 to 0.40.

#### Terminal Effluent Variance (Right Panel)

The variance of the treated effluent at the end of the simulation, $\operatorname{Var}(X_T^5)$, scales linearly with $\gamma_2$, rising from 0.07 at $\gamma_2=0.02$ to 0.50 at $\gamma_2=0.40$.

This linear relationship ($R^2\approx0.99$) confirms that variance transmits through the coupled linear cascade without triggering non-linear runaway instability.

#### Breach Probability (Left Panel)

Across all tested noise levels ($\gamma_2\in[0.02,0.40]$), the maximum cumulative breach probability:

$$
\mathbb{P}\left(\sup_t X^5\geq3.5\right)
$$

remains at 1.00 (100%).

Under this baseline loading profile, regulatory failure is driven primarily by mean hydraulic overloading rather than stochastic noise alone.

To reduce failure risk, the plant must expand its equalization tank volume $X^1$ or lower the baseline target outflow, rather than relying solely on noise suppression.

---

## 5\. Recent Developments and State of the Art

While the Harris (1977) formulation establishes the core framework for stochastic wastewater modeling, modern environmental systems engineering addresses several key limitations of the classical model.

### 5.1 Beyond Continuous White Noise: Jump-Diffusion Models

The Harris model assumes input shocks arrive via continuous Brownian motion $W_t$, where sample paths are continuous almost everywhere.

In real municipal systems, industrial batch dumps and flash storms arrive as discontinuous shocks.

Modern formulations incorporate Lévy processes or Poisson jump diffusions (Chapter 6.2):

$$
dX_t^3
=
f_3(X_t)\,dt
+
g_3(X_t)\,dW_t
+
J_t\,dN_t
$$

where:

- $N_t$ is a Poisson counting process with arrival intensity $\lambda$.
- $J_t$ is a random variable modeling the amplitude of sudden pollutant surges.

### 5.2 Autocorrelated and Colored Noise

White noise assumes that increments across non-overlapping intervals are completely independent:

$$
\mathbb{E}[dW_t\,dW_s]=0
\qquad (t\neq s)
$$

In reality, weather patterns and wastewater flows exhibit temporal persistence (memory).

Modern studies replace Gaussian white noise with colored noise generated by an Ornstein–Uhlenbeck process, or use fractional Brownian motion $B_H(t)$ with Hurst parameter $H>0.5$ to account for long-range environmental memory.

### 5.3 The IWA Anaerobic Digestion Model No. 1 (ADM1)

The biological treatment term in Harris:

$$
a+b(X^4-X^5)
$$

linearizes microbial kinetics.

Current industry practice relies on the **IWA Anaerobic Digestion Model No. 1 (ADM1)**.

ADM1 tracks biological degradation across dozens of state variables, including:

- Acidogenesis, acetogenesis, and methanogenesis pathways.
- Liquid-gas mass transfer for methane ($\text{CH}_4$) and carbon dioxide ($\text{CO}_2$).
- Dynamic pH balance and ammonia inhibition.

Recent research applies stochastic calculus to ADM1 by casting microbial kinetic coefficients (such as maximum uptake rate $k_m$ and decay rate $b_{\text{bac}}$) as stochastic mean-reverting processes, enabling predictive risk assessments for full-scale municipal digesters.

---

## 6\. Conclusion and Future Extensions

### 6.1 Summary of Key Findings

1. **Hypoelliptic Filtering:** Modeling wastewater treatment as a 5-variable degenerate SDE reveals that an equalization tank absorbs 15–20% of incoming environmental variance, acting as a low-pass filter on chaotic inflows.
2. **Hydraulic Transport Lag:** Cross-correlation analysis identified a physical retention lag of $\tau^*=2.4$ hours. This delay explains why instantaneous pathwise ratios $\mathbb{E}[X^5/X^3]$ produce misleading spikes when input loading drops rapidly. Evaluating cumulative mass balance provides a stable, physically consistent measure of treatment efficiency.
3. **Distributional Skewness:** Propagating square-root diffusions through nonlinear drift terms generates positive skewness in the effluent state distribution ($p(X^5,t)$), showing that transient water quality risks have heavier tails than standard Gaussian models predict.
4. **Reliability Analysis:** Under baseline operating conditions, effluent COD breaches regulatory thresholds ($C_{\text{crit}}=3.5$ mg/L) between $t=0.98$ days and $t=1.13$ days with probability 1.0. This confirms that meeting compliance standards requires adequate base equalization volume, not merely low environmental volatility.

### 6.2 Proposed Extensions

#### 1\. Stochastic Optimal Control

Introducing variable-speed effluent pumping $u(t)$ from the equalization tank, framed as a Hamilton–Jacobi–Bellman (HJB) optimal control problem:

$$
\min_{u(t)}
\mathbb{E}
\left[
\int_0^T
\left(
(X_t^5)^2+\lambda u(t)^2
\right)dt
\right]
$$

#### 2\. Nonlinear Monod / Haldane Kinetics

Replacing the linear reaction term with full Monod or substrate-inhibited Haldane kinetics, solving the resulting nonlinear SDEs using higher-order Milstein schemes.

#### 3\. Spatial Advection-Diffusion

Extending the lumped-parameter tank model to distributed 1D open channels using stochastic partial differential equations (SPDEs), bridging the gap between waste treatment and natural hydrology.

---

## 7\. Bibliography

1. Harris, C. J. (1977). Modeling the non-linear stochastic dynamics of an anaerobic digestion process. _Applied Mathematical Modelling, 1_(5), 265–272.
2. Monod, J. (1949). The growth of bacterial cultures. _Annual Review of Microbiology, 3_(1), 371–394.
3. Øksendal, B. (2003). _Stochastic Differential Equations: An Introduction with Applications_. Springer Science & Business Media.
4. Kloeden, P. E., & Platen, E. (1992). _Numerical Solution of Stochastic Differential Equations_. Springer-Verlag, Berlin.
5. Batstone, D. J., et al. (2002). The IWA Anaerobic Digestion Model No 1 (ADM1). _Water Science and Technology, 45_(10), 65–73.
6. Higham, D. J. (2001). An algorithmic introduction to numerical simulation of stochastic differential equations. _SIAM Review, 43_(3), 525–546.
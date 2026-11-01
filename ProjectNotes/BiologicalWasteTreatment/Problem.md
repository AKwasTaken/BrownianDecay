# The Core Problem
The goal is to use microorganisms (bacteria) to eat and break down urban wastewater (a process tracing back to the "[[BiologicalWaste#Monod model|Monod model]]" for bacterial growth).
The problem is that urban wastewater is highly unpredictable. Sudden rainfall, industrial dumping, or morning domestic routines cause massive spikes in both the amount of water and the toxicity of the sewage.
To handle this, Harris modeled an **anaerobic digester** using two connected tanks:
1. **The Equalization Tank:** A holding tank that receives the raw, chaotic sewage first. Its job is to mix and stabilize the flow so the bacteria don't get overwhelmed by a sudden toxic shock.
2. **The Reactor Tank (Continuously Stirred Tank Reactor):** The tank where the bacteria actually live and break down the pollutants.
![[Pasted image 20261101171316.png]]
__________
# The Five State Variables

To model this system mathematically, Harris tracks five specific variables (states) over time:
- **$X^1$**: The volume of the liquid in the Equalization Tank.
- **$X^2$**: The incoming flow rate of the raw sludge.
- **$X^3$**: The incoming concentration of the substrate (the raw pollutant/sewage).
- **$X^4$**: The concentration of the pollutant _leaving_ the Equalization Tank and entering the Reactor Tank.
- **$X^5$**: The concentration of the pollutant _leaving_ the Reactor Tank (this is the final, treated water).

| **Variable** | **What it measures** | **Standard Mathematical Units**         | **Example of a real-world value**            |
| ------------ | -------------------- | --------------------------------------- | -------------------------------------------- |
| **$X^1$**    | Volume in Eq. Tank   | Cubic meters ($m^3$) or Liters (L)      | $500 \text{ } m^3$ (Tank is half full)       |
| **$X^2$**    | Incoming Flow Rate   | Cubic meters per day ($m^3/\text{day}$) | $1,200 \text{ } m^3/\text{day}$ (Heavy flow) |
| **$X^3$**    | Incoming Dirtiness   | Milligrams per Liter (mg/L)*            | $800 \text{ mg/L}$ (Very strong sewage)      |
| **$X^4$**    | Mixed Dirtiness      | Milligrams per Liter (mg/L)             | $650 \text{ mg/L}$ (Smoothed out by tank)    |
| **$X^5$**    | Final Dirtiness      | Milligrams per Liter (mg/L)             | $50 \text{ mg/L}$ (Successfully treated)     |
*Engineers usually measure dirtiness using COD (Chemical Oxygen Demand), which tracks how much oxygen is needed to break down the waste.*

__________
# The Math: Breaking Down the Equations
Harris uses a 5-dimensional mathematical system to track these variables. Because the incoming sewage is unpredictable, the first two equations use **Itô Stochastic Differential Equations (SDEs)**, which include mathematical "noise."
#### Part A: The Random Inputs (Equations 6.1 & 6.2)
These two equations model the turbulent, chaotic nature of the raw sewage entering the plant.
$$(6.1) \qquad dX_t^2 = \left(\alpha_1 + 1 - \beta_1 \frac{X_t^2}{X_t^1}\right) dt + \sqrt{2\gamma_1 \frac{X_t^2}{X_t^1}} dW_t^1$$
$$(6.2) \qquad dX_t^3 = \frac{(\alpha_2 + 1 - \beta_2 X_t^3)}{X_t^1} dt + \sqrt{2\gamma_2 \frac{X_t^3}{X_t^1}} dW_t^2$$
- **What this means:** $dX_t^2$ tracks the changing flow rate, and $dX_t^3$ tracks the changing pollutant concentration.
- **The two halves:** In both equations, the part multiplied by $dt$ represents the predictable, average baseline behaviour. The part multiplied by $dW_t$ represents the unpredictable "white noise" or random spikes (like a sudden rainstorm).
- **The parameters:** $\alpha, \beta,$ and $\gamma$ are just positive numbers used to tune how severe the random fluctuations are.
#### Part B: The Equalization Tank (Equations 6.3 & 6.4)
Once the water is in the holding tank, we use standard **Ordinary Differential Equations (ODEs)** based on mass balance.
$$(6.3) \qquad dX_t^1 = \left(\frac{\alpha_1 + 1}{\beta_1} - X_t^2\right) dt$$
$$(6.4) \qquad dX_t^4 = (X_t^3 - X_t^4) \frac{X_t^2}{X_t^1} dt$$
- **What this means:** Equation 6.3 tracks how the volume of the holding tank goes up and down based on the incoming flow ($X^2$). Equation 6.4 tracks how the concentration of the pollutant gets diluted or concentrated before it leaves the holding tank ($X^4$).
- **No explicit noise:** Notice there is no $dW_t$ (noise) term here. The randomness is already accounted for inside the $X^2$ and $X^3$ variables from the previous step.
#### Part C: The Reactor Tank (Equation 6.5)
This is where the biological treatment happens.
$$(6.5) \qquad dX_t^5 = \left(a + b \left(X_t^4 - X_t^5\right)\right) dt$$
- **What this means:** This tracks the final pollutant concentration leaving the plant. It depends on what enters from the holding tank ($X^4$) and what is currently in the reactor ($X^5$).
- **The parameters:** $a$ and $b$ are biological factors. They represent how fast the micro-organisms grow and how efficiently they eat the waste (yield constant), as well as carrying over some noise parameters.
_________
# Conclusion
The 5-equation model functions as a nonlinear stochastic differential equation (SDE). Its most notable feature is how it isolates randomness: only the first two equations—which represent the system's inputs—contain explicit noise. Because the treatment tanks do not generate any new noise on their own and simply react to the incoming water, three of the five rows in the system's diffusion matrix are exactly zero.

To analyse how this system behaves, Harris ran a computer simulation of the equations. His goal was to evaluate the plant's overall performance under these chaotic conditions. To do this, he plotted the expected ratio of the final treated water ($X_t^5$) to the raw incoming water ($X_t^3$) ($E(X_t^5/X_t^3)$) over time. This graph provided a clear, measurable picture of how effectively the plant was reducing pollution.
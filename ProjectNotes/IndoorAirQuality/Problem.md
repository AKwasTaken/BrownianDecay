# The Core Problem & The Physical Setup
The goal is to figure out how the concentration of a pollutant specifically Carbon Dioxide ($CO_2$) changes inside a single room, and how to calibrate an air conditioning (HVAC) unit to filter it out.
In a perfect, laboratory-controlled world, you would know exactly how much $CO_2$ is in the room, how much is outside, and exactly how fast the AC pushes air around.
___
# The Baseline Math (The Perfect World)
Traditionally, engineers used a deterministic mass-balance equation to model this:

$$
(6.9) \qquad v \dot{C}^i = -k \left(q^0 + f^1 q^1 + q^2\right) C^i + S - R + k \left(q^0(1 - f^0) + q^2\right)C^0
$$

This looks intimidating, but it is just an accounting ledger for air. Let's break down the variables:
- **The Basics:**
    - $v \dot{C}^i$: The room's volume ($v$) multiplied by the changing indoor concentration of $CO_2$ ($\dot{C}^i$). This is what we are trying to solve for.
    - $C^0$: The outdoor $CO_2$ concentration.
- **The Airflow ($q$ variables):** Air moves in different ways.
    - $q^0$ is the "makeup air" (fresh air intentionally brought in by the AC).
    - $q^1$ is recirculated air (air from the room, pushed back through the AC).
    - $q^2$ is infiltration (air sneaking in through cracks in windows or under doors).
- **The Filters ($f$ variables):**
    - $f^0$ is how well the AC filters the outdoor makeup air.
    - $f^1$ is how well it filters the recirculated indoor air.
- **Sources, Sinks, and Mixing:**
    - $S$: The indoor source of $CO_2$ (e.g., people breathing).
    - $R$: The indoor sink (things removing $CO_2$, like an indoor plant, though usually negligible for $CO_2$).
    - $k$: A mixing factor. Air doesn't mix perfectly instantly; this accounts for that inefficiency.

|**Symbol**|**Description**|**Typical Units (SI)**|**Mathematical Form / Definition**|
|---|---|---|---|
|$v$|Room or zone volume|$m^3$|$L \times W \times H$|
|$C^i$|Indoor CO$_2$ concentration|$mg/m^3$ (or ppm)|$M_{CO_2} / v$ (Mass of CO$_2$ / Volume)|
|$\dot{C}^i$|Rate of change of indoor CO$_2$ concentration|$mg/m^3/hr$|$\frac{dC^i}{dt}$|
|$v \dot{C}^i$|Total mass accumulation rate of indoor CO$_2$|$mg/hr$|$v \frac{dC^i}{dt}$|
|$C^0$|Outdoor (ambient) CO$_2$ concentration|$mg/m^3$ (or ppm)|Assumed constant or time-series $C^0(t)$|
|$q^0$|Volumetric flow rate of makeup (fresh) air|$m^3/hr$|$\frac{dV_{fresh}}{dt}$|
|$q^1$|Volumetric flow rate of recirculated air|$m^3/hr$|$\frac{dV_{recirc}}{dt}$|
|$q^2$|Volumetric flow rate of infiltration air|$m^3/hr$|$\frac{dV_{infil}}{dt}$|
|$f^0$|Filter efficiency for makeup air|Dimensionless ($0 \le f \le 1$)|$1 - \frac{C_{after\_filter}}{C_{before\_filter}}$|
|$f^1$|Filter efficiency for recirculated air|Dimensionless ($0 \le f \le 1$)|$1 - \frac{C_{after\_filter}}{C_{before\_filter}}$|
|$S$|Indoor source generation rate (e.g., human respiration)|$mg/hr$|$\sum (\text{emission rate per person} \times \text{occupants})$|
|$R$|Indoor sink removal rate (non-ventilation losses)|$mg/hr$|$\frac{dM_{removed}}{dt}$|
|$k$|Mixing factor (ventilation effectiveness)|Dimensionless ($0 < k \le 1$)|$\frac{C_{exhaust} - C_{supply}}{C^i - C_{supply}}$|

> Note on units: In indoor air quality (IAQ) mass balance equations, if concentrations ($C^i$, $C^0$) are measured in parts-per-million (ppm), they are typically treated as dimensionless volumetric fractions ($m^3$ of CO$_2$ per $m^3$ of air). In that case, the source ($S$) and sink ($R$) rates must be expressed in volumetric terms ($m^3/hr$) rather than mass terms ($mg/hr$) to balance the equation.
___
# Injecting the Real-World Noise
In reality, a building is not a sealed lab. People enter and exit (changing the breathing rate), doors open and close (changing airflow), and a delivery truck idling outside changes the outdoor air quality.
Researchers Haghighat, Fazio, and Unny decided that four of these variables were too unpredictable to be static numbers. They replaced them with an average value (marked with a bar, like $\bar{S}$) plus a random noise fluctuation (marked with $\xi$):
1. **Indoor Source ($S$):** Becomes $\bar{S} + \xi^{(1)}$ (e.g., people unexpectedly entering the room).
2. **Outdoor Concentration ($C^0$):** Becomes $\bar{C}^0 + \xi^{(2)}$ (e.g., a truck driving by).
3. **Makeup Air Flow ($q^0$):** Becomes $\bar{q}^0 + \xi^{(3)}$ (e.g., mechanical fluctuations in the AC fan).
4. **Infiltration ($q^2$):** Becomes $\bar{q}^2 + \xi^{(4)}$ (e.g., wind gusts pushing air through window cracks).

$$
S_t = \bar{S} + \xi_t^{(1)}, \qquad C_t^0 = \bar{C}^0 + \xi_t^{(2)}
$$


$$
q_t^0 = \bar{q}^0 + \xi_t^{(3)}, \qquad q_t^2 = \bar{q}^2 + \xi_t^{(4)}
$$

# The SDE and Derivation
These stochastic parameters were substituted back into the baseline differential equation. During this derivation, the higher-order cross-product of the noise terms $\xi_t^{(2)}\xi_t^{(3)}$ (the interaction between outdoor concentration noise and make-up air flow noise) was neglected to yield a tractable Ito stochastic differential equation (SDE):

$$
v dC_t^i = g(C^i, t) dt + H(C^i, t) d\xi_t
$$

Here, $\xi_t = (\xi_t^{(1)}, \xi_t^{(2)}, \xi_t^{(3)}, \xi_t^{(4)})$ represents the combined noise processes (interpreted as independent Wiener processes).
The deterministic drift term, $g(C^i, t)$, is a scalar function based on the mean values:

$$
g(C^i, t) = -k \left(\bar{q}^0 + f^1 q^1 + \bar{q}^2\right) C^i + \bar{S} - R + k \left(\bar{q}^0(1 - f^0) + \bar{q}^2\right)\bar{C}^0
$$

The stochastic diffusion term is governed by $H(C^i, t)$, a row vector function mathematically capturing how each specific noise term shocks the system. Its transpose is explicitly defined as:

$$
H(C^i, t)^T = \begin{pmatrix} 1 \\ k \left(\bar{q}^0 (1 - f^0) + \bar{q}^2\right) \\ k(\bar{q}^0(1 - f^0) - C^i)\bar{C}^0 \\ k(\bar{C}^0 - C^i) \end{pmatrix}
$$

___
# The Final SDE (Equation 6.10)
They substituted these noisy variables back into the original Equation 6.9.
_(A quick math trick they used: when you multiply two random noise variables together, like $\xi^{(2)} \times \xi^{(3)}$, the resulting number is so unimaginably small that they just crossed it out to keep the math clean)._
The result is the final **Itô Stochastic Differential Equation**:

$$
(6.10) \qquad v dC_t^i = g(C^i, t) dt + H(C^i, t) d\xi_t
$$

This beautifully separates the predictable physics from the chaotic noise:
- **The Drift Term, $g(C^i, t) dt$:** This is just the original equation using the _average_ baseline values. It represents what the $CO_2$ level is mathematically expected to do over time.
- **The Diffusion Term, $H(C^i, t) d\xi_t$:** This is a vector matrix ($H$) multiplied by the four sources of noise ($\xi_t$). It mathematically calculates how much a sudden gust of wind, a passing truck, or an extra person breathing will instantly knock the $CO_2$ concentration off its expected path.
The derivation for the Itô SDE for air quality is given []
_____________
# The Bigger Picture
Just like the hydrologists chained their equations to model the Great Lakes, the authors noted that you can link these equations together to model air quality moving through **several different connected rooms** in a building.
Furthermore, the text ends by pointing out the versatility of this math. Another research team took this exact same mathematical structure—a deterministic baseline plus a matrix of noise—and used it to model **thermal variations** (temperature changes) in a building, rather than pollutants.
___
The Research Paper by F. HAGHIGHAT, P. FAZIO, T. E. UNNY on [[AirQualityPaper_Unny.pdf|A Predictive Stochastic Model For Indoor Air Quality]].
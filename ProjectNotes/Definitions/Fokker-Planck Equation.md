### The Fokker-Planck Equation (The Math)

The Fokker-Planck equation is a partial differential equation (PDE) that tracks how a probability distribution $p(x,t)$ changes over time. In a simplified 1D form, it looks like this:

$$\frac{\partial p(x,t)}{\partial t} = -\frac{\partial}{\partial x} \big[ D_1(x) p(x,t) \big] + \frac{\partial^2}{\partial x^2} \big[ D_2(x) p(x,t) \big]$$

- **$p(x,t)$**: The probability of the reservoir having exactly $x$ amount of water at time $t$.
    
- **$D_1(x)$ (Drift Term)**: The deterministic forces pushing the water level up or down (e.g., constant rain and evaporation).
    
- **$D_2(x)$ (Diffusion Term)**: The random, unpredictable fluctuations (noise) causing the water levels to vary.
    
- **$\frac{\partial p}{\partial t}$**: How the overall probability curve is shifting or flattening out as time ticks forward.
    

___
In the [[HydrologyPaper.pdf|paper]], the authors derive a specific Fokker-Planck equation to match the linear reservoir model (Equation 26):

$$\frac{\partial p}{\partial t} = -\frac{\partial}{\partial y} \big[ p(K_1y + K_2) \big] + \frac{1}{2}\frac{\partial^2}{\partial y^2} \big[ p(K_3y^2 + K_4y + K_5) \big]$$

Here is how it connects directly to the project variables:

$y$ (Outflow): The authors use $y$ instead of $x$ to represent the reservoir outflow.

$K_1y + K_2$ (The Drift): This encapsulates the deterministic physics—the baseline average precipitation, evaporation, and human abstractions.

$K_3y^2 + K_4y + K_5$ (The Diffusion): This represents the combined variance (noise) from the fluctuating weather. The $K$ coefficients are built directly from the $\sigma_1, \sigma_2,$ and $\sigma_3$ noise values.

$\frac{\partial p}{\partial t} = 0$ (Steady-State): If you wait long enough, the reservoir system reaches equilibrium. The probability curve stops shifting, meaning the rate of change over time ($\frac{\partial p}{\partial t}$) drops to exactly zero.
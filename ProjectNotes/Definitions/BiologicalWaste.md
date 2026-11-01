# Monod model
The **Monod model** is a mathematical model used in microbiology and biochemical engineering to describe the growth rate of microorganisms. Proposed by Jacques Monod in 1949, it establishes how the specific growth rate of a microbial population depends on the concentration of a limiting nutrient (substrate) in their environment.

It is the foundational equation for designing bioreactors, modelling wastewater treatment processes, and understanding environmental microbial ecology.

## The Equation

The model is expressed by a hyperbolic equation that looks very similar to the Michaelis-Menten equation used in enzyme kinetics:


$$
\mu = \mu_{max} \frac{S}{K_s + S}
$$


Where:

- **$\mu$ (mu):** The specific growth rate of the microorganisms (usually measured in inverse time, e.g., $h^{-1}$). This is how fast the population is multiplying.
    
- **$\mu_{max}$:** The maximum specific growth rate. This is the theoretical top speed at which the cells can reproduce when food is completely unlimited.
    
- **$S$:** The concentration of the limiting substrate (the essential nutrient that is in the shortest supply, like glucose or nitrogen).
    
- **$K_s$:** The half-velocity constant (or half-saturation constant). This is the exact substrate concentration at which the growth rate is exactly half of $\mu_{max}$ ($\mu = \frac{\mu_{max}}{2}$).
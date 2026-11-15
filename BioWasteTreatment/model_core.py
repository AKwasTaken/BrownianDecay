import os
import numpy as np

# Output directory matching LaTeX \includegraphics{img/...}
OUTPUT_DIR = "img"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Time grid & Monte Carlo settings
T = 5.0
N_steps = 2000
dt = T / N_steps
sq_dt = np.sqrt(dt)
t = np.linspace(0.0, T, N_steps + 1)
N_paths = 3000

# Model parameters
alpha1, beta1, gamma1 = 1.0, 1.2, 0.15
alpha2, beta2, gamma2 = 1.0, 1.0, 0.20
D = 1.35
k = 0.25
a = 0.05
b = D + k  # 1.60

# Initial conditions [X1_0, X2_0, X3_0, X4_0, X5_0]
X1_0 = 10.0
X2_0 = ((alpha1 + 1.0) / beta1) * X1_0  # Equilibrium inflow: 16.667
X3_0 = 5.0
X4_0 = 4.0
X5_0 = 1.0

# Regulatory & Controller thresholds
C_crit = 3.5
C_target = 2.6
k_c = 6.0
u_min = 0.05


def simulate_passive(seed=42):
    """Projected Euler-Maruyama solver for the passive 5D Harris system."""
    rng = np.random.default_rng(seed)
    X = np.zeros((5, N_steps + 1, N_paths))
    X[:, 0, :] = np.array([X1_0, X2_0, X3_0, X4_0, X5_0])[:, None]

    for n in range(N_steps):
        x1 = np.maximum(X[0, n, :], 1e-4)
        x2 = np.maximum(X[1, n, :], 0.0)
        x3 = np.maximum(X[2, n, :], 0.0)
        x4 = np.maximum(X[3, n, :], 0.0)
        x5 = np.maximum(X[4, n, :], 0.0)

        Z1 = rng.standard_normal(N_paths)
        Z2 = rng.standard_normal(N_paths)

        diff2 = np.sqrt(2.0 * gamma1 * x2 / x1) * sq_dt * Z1
        diff3 = np.sqrt(2.0 * gamma2 * x3 / x1) * sq_dt * Z2

        q_out = ((alpha1 + 1.0) / beta1) * x1

        X[0, n + 1, :] = np.maximum(x1 + (x2 - q_out) * dt, 1e-4)
        X[1, n + 1, :] = np.maximum(x2 + (alpha1 + 1.0 - beta1 * (x2 / x1)) * dt + diff2, 0.0)
        X[2, n + 1, :] = np.maximum(x3 + ((alpha2 + 1.0 - beta2 * x3) / x1) * dt + diff3, 0.0)
        X[3, n + 1, :] = np.maximum(x4 + ((x3 - x4) * (x2 / x1)) * dt, 0.0)
        X[4, n + 1, :] = np.maximum(x5 + (D * x4 - b * x5 - a) * dt, 0.0)

    return X


def simulate_controlled(seed=42):
    """Projected Euler-Maruyama solver with coupled sigmoidal throttle valve u(X5)."""
    rng = np.random.default_rng(seed)
    X = np.zeros((5, N_steps + 1, N_paths))
    X[:, 0, :] = np.array([X1_0, X2_0, X3_0, X4_0, X5_0])[:, None]

    for n in range(N_steps):
        x1 = np.maximum(X[0, n, :], 1e-4)
        x2 = np.maximum(X[1, n, :], 0.0)
        x3 = np.maximum(X[2, n, :], 0.0)
        x4 = np.maximum(X[3, n, :], 0.0)
        x5 = np.maximum(X[4, n, :], 0.0)

        Z1 = rng.standard_normal(N_paths)
        Z2 = rng.standard_normal(N_paths)

        diff2 = np.sqrt(2.0 * gamma1 * x2 / x1) * sq_dt * Z1
        diff3 = np.sqrt(2.0 * gamma2 * x3 / x1) * sq_dt * Z2

        u = u_min + (1.0 - u_min) / (1.0 + np.exp(k_c * (x5 - C_target)))
        q_out = u * ((alpha1 + 1.0) / beta1) * x1

        X[0, n + 1, :] = np.maximum(x1 + (x2 - q_out) * dt, 1e-4)
        X[1, n + 1, :] = np.maximum(x2 + (alpha1 + 1.0 - beta1 * (x2 / x1)) * dt + diff2, 0.0)
        X[2, n + 1, :] = np.maximum(x3 + ((alpha2 + 1.0 - beta2 * x3) / x1) * dt + diff3, 0.0)
        X[3, n + 1, :] = np.maximum(x4 + ((x3 - x4) * (x2 / x1)) * dt, 0.0)
        X[4, n + 1, :] = np.maximum(x5 + (D * u * x4 - (D * u + k) * x5 - a) * dt, 0.0)

    return X

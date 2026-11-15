import os
import numpy as np
import matplotlib.pyplot as plt
from model_core import OUTPUT_DIR, N_steps, dt, t, simulate_passive

X = simulate_passive()

# 1. Ensemble Variance Transmission Ratio Var(X4) / Var(X3)
var_x3 = np.var(X[2, :, :], axis=1)
var_x4 = np.var(X[3, :, :], axis=1)
damping_ratio = np.zeros_like(var_x3)
damping_ratio[1:] = var_x4[1:] / np.maximum(var_x3[1:], 1e-8)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, damping_ratio, color="crimson", lw=2)
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"$\rho_{\mathrm{damp}}(t)$")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "variance_damping_ratio.pdf"), format="pdf")
plt.close(fig)

# 2. Centered Cross-Correlation Lag R_34(tau)
x3_c = X[2, :, :] - np.mean(X[2, :, :], axis=0, keepdims=True)
x4_c = X[3, :, :] - np.mean(X[3, :, :], axis=0, keepdims=True)

max_lag = int(N_steps / 2)
lags = np.arange(-max_lag, max_lag + 1)
ccorr = np.zeros(len(lags))

for idx, lag in enumerate(lags):
    if lag < 0:
        a_s, b_s = x3_c[-lag:, :], x4_c[:lag, :]
    elif lag > 0:
        a_s, b_s = x3_c[:-lag, :], x4_c[lag:, :]
    else:
        a_s, b_s = x3_c, x4_c
    ccorr[idx] = np.mean(a_s * b_s) / (np.std(a_s) * np.std(b_s) + 1e-8)

lag_times = lags * dt
best_idx = np.argmax(ccorr)
tau_star = lag_times[best_idx]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(lag_times, ccorr, color="teal", lw=2)
ax.axvline(tau_star, color="black", linestyle=":", label=rf"Max Delay $\tau^* = {tau_star:.2f}$ days")
ax.set_xlabel(r"Lag $\tau$ (days)")
ax.set_ylabel("Normalized Correlation")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "cross_correlation_lag.pdf"), format="pdf")
plt.close(fig)

print(f"[3_variance_and_lag] Peak damping ratio = {np.max(damping_ratio):.2f}, terminal = {damping_ratio[-1]:.2f}")
print(f"[3_variance_and_lag] Optimal lag tau* = {tau_star:.2f} days ({tau_star * 24:.1f} hrs), corr = {ccorr[best_idx]:.2f}")

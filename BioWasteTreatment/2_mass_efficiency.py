import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from model_core import OUTPUT_DIR, t, alpha1, beta1, simulate_passive

X = simulate_passive()

q_out = ((alpha1 + 1.0) / beta1) * X[0, :, :]
mass_flow_in = X[1, :, :] * X[2, :, :]
mass_flow_out = q_out * X[4, :, :]

cum_mass_in = cumulative_trapezoid(mass_flow_in, t, axis=0, initial=0.0)
cum_mass_out = cumulative_trapezoid(mass_flow_out, t, axis=0, initial=0.0)

with np.errstate(divide="ignore", invalid="ignore"):
    eta_mass = (1.0 - (cum_mass_out / cum_mass_in)) * 100.0
eta_mass[0, :] = np.nan  # Undefined (0/0) at t = 0

mean_eta = np.nanmean(eta_mass, axis=1)

# 1. Cumulative Mass Removal Efficiency
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, mean_eta, color="darkgreen", lw=2.2)
ax.set_xlabel("Time (days)")
ax.set_ylabel("Waste Removed (%)")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "cumulative_mass_removal_efficiency.pdf"), format="pdf")
plt.close(fig)

# 2. Macro-Concentration Drop Ratio E[X5] / E[X3]
ratio_of_exp = np.mean(X[4, :, :], axis=1) / np.maximum(np.mean(X[2, :, :], axis=1), 1e-4)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, ratio_of_exp, color="indigo", lw=2.2, label=r"$\mathbb{E}[X^5_t] / \mathbb{E}[X^3_t]$")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Mean Ratio")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "macro_concentration_drop_ratio.pdf"), format="pdf")
plt.close(fig)

print(f"[2_mass_efficiency] Terminal cumulative removal = {mean_eta[-1]:.1f}%")
print(f"[2_mass_efficiency] Macro ratio range: {ratio_of_exp[0]:.2f} -> {np.max(ratio_of_exp):.2f}")

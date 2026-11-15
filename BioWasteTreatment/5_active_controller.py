import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
from model_core import (
    OUTPUT_DIR, t, alpha1, beta1, C_crit, C_target, k_c, u_min,
    simulate_passive, simulate_controlled
)

# Run both with identical seed=42 so noise paths match 1-to-1
X_base = simulate_passive(seed=42)
X_ctrl = simulate_controlled(seed=42)

# 1. Survival Comparison
surv_base = 1.0 - np.mean(np.maximum.accumulate(X_base[4, :, :], axis=0) >= C_crit, axis=1)
surv_ctrl = 1.0 - np.mean(np.maximum.accumulate(X_ctrl[4, :, :], axis=0) >= C_crit, axis=1)

# 2. Cumulative Mass Removal Efficiency (using true throttled outflow)
u_mat = u_min + (1.0 - u_min) / (1.0 + np.exp(k_c * (X_ctrl[4, :, :] - C_target)))
q_out_base = ((alpha1 + 1.0) / beta1) * X_base[0, :, :]
q_out_ctrl = u_mat * ((alpha1 + 1.0) / beta1) * X_ctrl[0, :, :]

mass_in_base = cumulative_trapezoid(X_base[1] * X_base[2], t, axis=0, initial=0.0)
mass_out_base = cumulative_trapezoid(q_out_base * X_base[4], t, axis=0, initial=0.0)
with np.errstate(divide="ignore", invalid="ignore"):
    eta_base = (1.0 - mass_out_base / mass_in_base) * 100.0
eta_base[0, :] = np.nan

mass_in_ctrl = cumulative_trapezoid(X_ctrl[1] * X_ctrl[2], t, axis=0, initial=0.0)
mass_out_ctrl = cumulative_trapezoid(q_out_ctrl * X_ctrl[4], t, axis=0, initial=0.0)
with np.errstate(divide="ignore", invalid="ignore"):
    eta_ctrl = (1.0 - mass_out_ctrl / mass_in_ctrl) * 100.0
eta_ctrl[0, :] = np.nan

# Plot 1: Effluent Concentration Comparison
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, np.mean(X_base[4, :, :], axis=1), color="crimson", lw=2.2, label="Passive Baseline")
ax.plot(t, np.mean(X_ctrl[4, :, :], axis=1), color="#006699", lw=2.2, label="Active Retention Control")
ax.axhline(C_crit, color="black", linestyle=":", lw=1.5, label=rf"Regulatory Limit $C_{{crit}} = {C_crit}$")
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"Effluent COD $X^5$ (mg/L)")
ax.legend(loc="upper left")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "control_effluent_comparison.pdf"), format="pdf")
plt.close(fig)

# Plot 2: Regulatory Survival Comparison
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, surv_base, color="crimson", lw=2.2, label="Passive Baseline")
ax.plot(t, surv_ctrl, color="#006699", lw=2.2, label="Active Retention Control")
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"Survival Probability $S(t) = P(\tau > t)$")
ax.legend(loc="lower left")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "control_survival_comparison.pdf"), format="pdf")
plt.close(fig)

# Plot 3: Cumulative Mass Removal Efficiency Comparison
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, np.nanmean(eta_base, axis=1), color="crimson", lw=2.2, label="Passive Baseline")
ax.plot(t, np.nanmean(eta_ctrl, axis=1), color="#006699", lw=2.2, label="Active Retention Control")
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"Cumulative Waste Removed $\eta_{\mathrm{mass}}$ (%)")
ax.set_ylim(-5, 100)
ax.legend(loc="upper right")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "control_efficiency_comparison.pdf"), format="pdf")
plt.close(fig)

print(f"[5_active_controller] 5-day Survival: Passive = {surv_base[-1]*100:.1f}%, Controlled = {surv_ctrl[-1]*100:.1f}%")
print(f"[5_active_controller] Terminal Removal: Passive = {np.nanmean(eta_base[-1]):.1f}%, Controlled = {np.nanmean(eta_ctrl[-1]):.1f}%")
print(f"[5_active_controller] Max Tank 1 Volume E[X1]: Passive = {np.max(np.mean(X_base[0], axis=1)):.1f}, Controlled = {np.max(np.mean(X_ctrl[0], axis=1)):.1f}")

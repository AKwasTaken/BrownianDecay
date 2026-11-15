import os
import numpy as np
import matplotlib.pyplot as plt
from model_core import OUTPUT_DIR, t, simulate_passive

X = simulate_passive()

# 1. Input Dynamics
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, X[1, :, :3], alpha=0.6, lw=1.2)
ax.plot(t, X[2, :, :3], alpha=0.6, lw=1.2, linestyle="--")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Flow / Concentration")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "input_dynamics.pdf"), format="pdf")
plt.close(fig)

# 2. Equalization Effect
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, X[2, :, 0], label="Raw Input ($X^3$)", color="crimson", alpha=0.7)
ax.plot(t, X[3, :, 0], label="Equalized ($X^4$)", color="navy", lw=2)
ax.set_xlabel("Time (days)")
ax.set_ylabel("Concentration (mg/L)")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "equalization_effect.pdf"), format="pdf")
plt.close(fig)

# 3. Effluent Quality
mean_x5 = np.mean(X[4, :, :], axis=1)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, X[4, :, :10], color="gray", alpha=0.3, lw=0.8)
ax.plot(t, mean_x5, color="darkgreen", lw=2.2, label=r"Ensemble Mean $\mathbb{E}[X^5]$")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Effluent Concentration ($X^5$) [mg/L]")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "effluent_quality.pdf"), format="pdf")
plt.close(fig)

# 4. Reduction Ratio E[X5 / X3]
ratio = X[4, :, :] / np.maximum(X[2, :, :], 1e-4)
mean_ratio = np.mean(ratio, axis=1)
pct_10 = np.percentile(ratio, 10, axis=1)
pct_90 = np.percentile(ratio, 90, axis=1)

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, mean_ratio, color="purple", lw=2.2, label=r"$\mathbb{E}[X^5_t / X^3_t]$")
ax.fill_between(t, pct_10, pct_90, color="purple", alpha=0.2, label="10th–90th Percentile")
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"Reduction Ratio $\mathbb{E}[X^5 / X^3]$")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "reduction_ratio.pdf"), format="pdf")
plt.close(fig)

print(f"[1_basic_dynamics] Peak E[X5] = {np.max(mean_x5):.2f} mg/L at t = {t[np.argmax(mean_x5)]:.2f} d")
print(f"[1_basic_dynamics] Terminal E[X5/X3] = {mean_ratio[-1]:.2f}")

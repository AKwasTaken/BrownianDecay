import os
import numpy as np
import matplotlib.pyplot as plt
from model_core import OUTPUT_DIR, t, C_crit, simulate_passive

X = simulate_passive()

breach_mask = X[4, :, :] >= C_crit
has_breached = np.any(breach_mask, axis=0)
first_breach_idx = np.argmax(breach_mask, axis=0)
fht_values = t[first_breach_idx[has_breached]]

running_max_x5 = np.maximum.accumulate(X[4, :, :], axis=0)
p_breach = np.mean(running_max_x5 >= C_crit, axis=1)
survival = 1.0 - p_breach

# 1. First Hitting Time Density
fig, ax = plt.subplots(figsize=(7, 4.5))
if len(fht_values) > 0:
    ax.hist(fht_values, bins=30, density=True, color="darkorange", alpha=0.7, edgecolor="black")
ax.set_xlabel(r"First Breach Time $\tau$ (days)")
ax.set_ylabel("Empirical Density")
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "first_hitting_time_density.pdf"), format="pdf")
plt.close(fig)

# 2. Survival & Failure Risk
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(t, survival, color="darkblue", lw=2, label=r"Survival $S(t) = P(\tau > t)$")
ax.plot(t, p_breach, color="firebrick", lw=2, linestyle="--", label=r"Failure $P(\tau \leq t)$")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Probability")
ax.legend()
ax.grid(True, linestyle="--", alpha=0.5)
fig.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "reliability_and_violation_risk.pdf"), format="pdf")
plt.close(fig)

print(f"[4_fht_and_reliability] Mean FHT = {np.mean(fht_values):.2f} days (median = {np.median(fht_values):.2f} d)")
print(f"[4_fht_and_reliability] 5-day Survival = {survival[-1]*100:.1f}%, Breach = {p_breach[-1]*100:.1f}%")

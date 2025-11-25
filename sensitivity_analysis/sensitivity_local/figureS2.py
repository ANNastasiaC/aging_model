results = {}

for name, bounds in param_bounds.items():
    low_val,  high_val = bounds

    low_key = name + "low"
    high_key = name + "high"

    low = sen_ratios[low_key]
    base = sen_ratios["baseline"]

    if name == 'k_c' or name =='k7':
        high = "Nan"
    else:
        high = sen_ratios[high_key]

    results[name] = (low, base, high)

# -------------------------
# Figure 2: Bar chart of sensitivities (with LaTeX labels)
# -------------------------

import matplotlib.pyplot as plt
import numpy as np

# --- LaTeX names matching dictionary (same as Fig. 1) ---
pretty_names = {
    "k_div_s": r"$k_{\mathrm{div},s}$",
    "k_div_d": r"$k_{\mathrm{div},d}$",
    "p_dif": r"$p_{\mathrm{dif}}$",
    "c_ros": r"$c_{\mathrm{ROS}}$",
    "f_e(t)": r"$f_{e}(t)$",
    "k2": r"$k_{\mathrm{dm}}$",
    "k_6": r"$p_{m,\mathrm{sen}}$",
    "k7": r"$k_{1}$",
    "k8": r"$k_{2}$",
    "k0": r"$k_{0}$",
    "k_a_s": r"$k_{a,s}$",
    "k_a_d": r"$k_{a,d}$",
    "k_a_sn": r"$k_{a,sn}$",
    "k_a_dn": r"$k_{a,dn}$",
    "r_rep_s": r"$r_{\mathrm{rep},s}$",
    "r_rep_d": r"$r_{\mathrm{rep},d}$",
    "k_c": r"$k_{c}$",
    "resistance_to_apoptosis": r"$k_{\mathrm{resa}}$",
    "r_immuno_surv": r"$r_{\mathrm{im}}$",
    "rs_ub": r"$b_{\mathrm{rs}}$",
    "k_q": r"$k_{q}$",
    "r_ra": r"$r_{\mathrm{rea}}$",
    "k_sasp": r"$k_{\mathrm{sasp}}$",
    "stem_cell_init": r"$f_{0,\mathrm{stem}}$"
}

# -------------------------
# Build plotting lists
# -------------------------

param_names = list(results.keys())
pretty_labels = [pretty_names.get(p, p) for p in param_names]
y_pos = np.arange(len(param_names))

bar_width = 0.32
baseline_val = sen_ratios["baseline"]

# -------------------------
# Plot
# -------------------------

fig, ax = plt.subplots(figsize=(8, 11), dpi=500)

for i, param in enumerate(param_names):
    low, base, high = results[param]

    # --- Low bar ---
    ax.barh(y_pos[i] - bar_width/2,
            low - baseline_val,
            height=bar_width,
            color='red',
            alpha=0.7,
            label="Low (0.1×)" if i == 0 else "")

    # --- High bar (skip exploding params) ---
    if not (param == 'k_c' or param == 'k7'):
        ax.barh(y_pos[i] + bar_width/2,
                high - baseline_val,
                height=bar_width,
                color='blue',
                alpha=0.7,
                label="High (10×)" if i == 0 else "")

# -------------------------
# Axis styling
# -------------------------

ax.set_yticks(y_pos)
ax.set_yticklabels(pretty_labels)

ax.axvline(0, color="black", linewidth=0.8)

ax.set_xlabel(r"Change in senescent fraction at age 80 (relative to baseline)")
# ax.set_title("Figure 2. Parameter sensitivity (±10× perturbation)")

ax.legend(loc="upper right")
plt.tight_layout()
plt.show()

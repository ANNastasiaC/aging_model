import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# ----------------------------------------
# A. MATCH PARAMETER NAMES IN CODE WITH THOSE USED IN PAPER
# ----------------------------------------
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

# ----------------------------------------
# B. STYLE SETTINGS
# ----------------------------------------
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=plt.cm.tab20.colors)

style_map = {
    "baseline": ("-", 3.0, 1.0),
    "low": ("--", 1.6, 0.7),
    "high": (":", 1.6, 0.7),
}

# ----------------------------------------
# C. PLOTTING
# ----------------------------------------

plt.figure(figsize=(10, 9), dpi = 500)

time = np.arange(len(sen_curves['baseline'])) * dt

# --- Baseline ---
ls, lw, alpha = style_map["baseline"]
plt.plot(time,
         sen_curves["baseline"],
         linestyle=ls, linewidth=lw, alpha=alpha,
         label="Baseline")

# --- Parameter variations ---
for param in param_bounds.keys():
    for variation in ["low", "high"]:
        key = f"{param}{variation}"
        if key not in sen_curves:
            continue

        ls, lw, alpha = style_map[variation]
        pretty_name = pretty_names.get(param, param)

        plt.plot(
            time,
            sen_curves[key],
            linestyle=ls, linewidth=lw, alpha=alpha,
            label=f"{pretty_name} ({variation})"
        )

# ----------------------------------------
# D. LABELS / LEGEND
# ----------------------------------------

plt.xlabel("Age (years)")
plt.ylabel("Senescent fraction")
#plt.title("Figure 1. Senescence ratio under parameter perturbations")

plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left', ncol=2)
plt.tight_layout()
plt.show()

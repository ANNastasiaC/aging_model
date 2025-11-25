import matplotlib.pyplot as plt
import numpy as np

# --- Pretty LaTeX names dictionary ---
pretty_names = {
    "k_div_s": r"$k_{\mathrm{div},s}$",
    "k_div_d": r"$k_{\mathrm{div},d}$",
    "p_dif": r"$p_{\mathrm{dif}}$",
    "c_ros": r"$c_{\mathrm{ROS}}$",
    "f_e(t)": r"$f_{e}hiii",
    "k2": r"$k_{\mathrm{dm}}$",
    "k_6": r"$p_{m,\mathrm{sen}}$",
    "k7": r"$k_1$",
    "k8": r"$k_2$",
    "k0": r"$k_0$",
    "k_a_s": r"$k_{a,s}$",
    "k_a_d": r"$k_{a,d}$",
    "k_a_sn": r"$k_{a,sn}$",
    "k_a_dn": r"$k_{a,dn}$",
    "r_rep_s": r"$r_{\mathrm{rep},s}$",
    "r_rep_d": r"$r_{\mathrm{rep},d}$",
    "k_c": r"$k_c$",
    "resistance_to_apoptosis": r"$k_{\mathrm{resa}}$",
    "r_immuno_surv": r"$r_{\mathrm{im}}$",
    "rs_ub": r"$b_{\mathrm{rs}}$",
    "k_q": r"$k_q$",
    "r_ra": r"$r_{\mathrm{rea}}$",
    "k_sasp": r"$k_{\mathrm{sasp}}$",
    "stem_cell_init": r"$f_{0,\mathrm{stem}}$"
}

# --- Sort PRCC by absolute value ---
sorted_idx = np.argsort(np.abs(list(PRCC.values())))[::-1]
sorted_names = np.array(list(PRCC.keys()))[sorted_idx]
sorted_values = np.array(list(PRCC.values()))[sorted_idx]

# --- Map to pretty names ---
pretty_labels = [pretty_names.get(name, name) for name in sorted_names]

# --- Plot ---
plt.figure(figsize=(6, 8), dpi = 500)  # reduced height
plt.barh(pretty_labels, sorted_values, color='skyblue', edgecolor='black', alpha=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.xlabel("Partial Rank Correlation Coefficient (PRCC)")
#plt.title("Sensitivity of 80-year Senescent Fraction")
plt.gca().invert_yaxis()  # largest PRCC on top
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.show()

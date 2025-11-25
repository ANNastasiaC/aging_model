# numba_simulation_full.py
import numpy as np
from numba import njit
import time

# -------------------------
# Model parameters
# -------------------------
initial = 1e4

# Time parameters
t_max = 5000.0
dt = 1.0
time_steps = int(t_max / dt) + 1

# Age parameters
a_max = 50.0
da = 1.0
age_steps = int(a_max / da) + 1

# Generation parameters
g_max = 80
dg = 1
gen_steps = int(g_max / dg) + 1


# Grids
t_grid = np.arange(0.0, t_max + dt, dt, dtype=np.float32)
a_grid = np.arange(0.0, a_max + da, da, dtype=np.float32)
g_grid = np.arange(0, g_max + dg, dg, dtype=np.int32)

# -------------------------
# Numba-compiled simulation
# -------------------------
@njit
def run_simulation_numba(initial, time_steps, age_steps, gen_steps, dt, k_div_s, k_div_d, p_dif, c_ros, f_e, k2, k_6, k7, k8, k0, k_a_s, k_a_d, k_a_sn, k_a_dn, r_rep_s, r_rep_d, k_c, resistance_to_apoptosis, r_immuno_surv, rs_ub, k_q, r_ra, k_sasp, stem_cell_init):

    # allocate arrays
    S = np.zeros((time_steps, age_steps, gen_steps), dtype=np.float32)
    D = np.zeros((time_steps, age_steps, gen_steps), dtype=np.float32)
    SN = np.zeros((time_steps, age_steps), dtype=np.float32)
    DN = np.zeros((time_steps, age_steps), dtype=np.float32)
    SQ = np.zeros((time_steps, age_steps, gen_steps), dtype=np.float32)

    # bookkeeping arrays
    total_d_sen_ms_array = np.zeros(time_steps, dtype=np.float32)
    total_d_sen_rs_array = np.zeros(time_steps, dtype=np.float32)
    total_d_sen_sasp_array = np.zeros(time_steps, dtype=np.float32)

    # initial conditions
    S[0, 0, 0] = initial * stem_cell_init
    D[0, 0, 0] = initial * (1-stem_cell_init)

    accumulated_dna_health = 1.0

    # time loop
    for t in range(1, time_steps):
        # compute r_dm for t
        r_dm_t = max(0.0, f_e + k2 * c_ros)
        accumulated_dna_health *= (1.0 - r_dm_t)

        # totals from previous timestep
        previous_sum = 0.0
        previous_sen = 0.0
        for ai in range(age_steps):
            for gi in range(gen_steps):
                previous_sum += S[t-1, ai, gi] + D[t-1, ai, gi] + SQ[t-1, ai, gi]
            previous_sum += SN[t-1, ai] + DN[t-1, ai]
            previous_sen += SN[t-1, ai] + DN[t-1, ai]

        if previous_sum <= 0.0:
            previous_sum = 1.0

        sen_ratio = previous_sen / previous_sum

        # crowding factor f_c using inlined upper_limit
        year = (t * dt) / 50.0
        weight = 70.0 - 80.0 / (1.0 + np.exp((year - 7.0) / 5.0))
        Nb = initial / 5.0 * weight
        f_c = Nb / previous_sum

        # r_dm t-1 for senescence computations
        r_dm_prev = max(0.0, f_e + k2 * c_ros)
        unrepaired_s = r_dm_prev * (1.0 - r_rep_s)
        unrepaired_d = r_dm_prev * (1.0 - r_rep_d)

        r_ms_s_val = 1.0 - (1.0 - k_6) ** unrepaired_s
        r_ms_d_val = 1.0 - (1.0 - k_6) ** unrepaired_d

        # update quiescent pool
        for ai in range(age_steps):
            for gi in range(gen_steps):
                SQ[t, ai, gi] = SQ[t-1, ai, gi] * (1.0 - r_ms_s_val * dt)

        # update senescent age arrays
        for a in range(1, age_steps):
            a_by_damage_sn = 1.0 - (1.0 - k_a_sn) ** unrepaired_s
            a_by_damage_dn = 1.0 - (1.0 - k_a_dn) ** unrepaired_d

            r_a_sn_val = max(0.0, (a_by_damage_sn + k_c * (1.0 - f_c)) * resistance_to_apoptosis + r_immuno_surv * accumulated_dna_health)
            r_a_dn_val = max(0.0, (a_by_damage_dn + k_c * (1.0 - f_c)) * resistance_to_apoptosis + r_immuno_surv * accumulated_dna_health)

            SN[t, a] = SN[t-1, a-1] - SN[t-1, a-1] * r_a_sn_val * dt
            DN[t, a] = DN[t-1, a-1] - DN[t-1, a-1] * r_a_dn_val * dt

        # reset bookkeeping counters for this timestep
        total_d_sen_ms = 0.0
        total_d_sen_rs = 0.0
        total_d_sen_sasp = 0.0

        # loop over ages and generations
        for a in range(1, age_steps):
            for g in range(gen_steps):
                # ----- STEM CELLS -----
                a_by_damage_s = 1.0 - (1.0 - k_a_s) ** unrepaired_s
                r_a_s_val = max(0.0, a_by_damage_s + k_c * (1.0 - f_c))
                r_s_sasp_val = sen_ratio * k_sasp
                total_rate_s = r_ms_s_val + r_a_s_val + r_s_sasp_val
                r_div_s_val = k_div_s * f_c
                r_q_val = max(0.0, k_q * (1.0 - f_c))
                r_reactivate_val = r_ra

                prev_S = S[t-1, a-1, g]
                S_val = prev_S * (1.0 - total_rate_s * dt) * (1.0 - r_div_s_val * dt - r_q_val * dt) + r_reactivate_val * SQ[t-1, a, g] * dt
                S[t, a, g] = S_val

                SN[t, a] += (r_ms_s_val + r_s_sasp_val) * prev_S * dt + r_ms_s_val * SQ[t-1, a-1, g]
                total_d_sen_sasp += r_s_sasp_val * D[t-1, a-1, g] * dt
                total_d_sen_ms += r_ms_d_val * D[t-1, a-1, g]

                # update quiescent
                SQ[t, a, g] += prev_S * (1.0 - total_rate_s * dt) * r_q_val * dt - r_reactivate_val * SQ[t-1, a, g] * dt

                # dividing stem cells and replicative senescence
                dividing_s_cells = prev_S * (1.0 - total_rate_s * dt) * r_div_s_val
                k1s = k7 + k8 * unrepaired_s
                r_rs_s_val = rs_ub / (1.0 + np.exp(-k1s * g + k0))
                if g + 1 < gen_steps:
                    S[t, 0, g+1] += dividing_s_cells * 2.0 * (1.0 - p_dif) * (1.0 - r_rs_s_val)
                D[t, 0, 0] += dividing_s_cells * 2.0 * p_dif * (1.0 - r_rs_s_val)
                SN[t, 0] += dividing_s_cells * 2.0 * (1.0 - p_dif) * r_rs_s_val
                DN[t, 0] += dividing_s_cells * 2.0 * p_dif * r_rs_s_val
                total_d_sen_rs += dividing_s_cells * 2.0 * p_dif * r_rs_s_val

                # ----- DIFFERENTIATED CELLS -----
                r_a_d_val = max(0.0, (1.0 - (1.0 - k_a_d) ** unrepaired_d) + k_c * (1.0 - f_c))
                total_rate_d = r_ms_d_val + r_s_sasp_val + r_a_d_val
                prev_D = D[t-1, a-1, g]
                D_val = prev_D * (1.0 - total_rate_d * dt) * (1.0 - k_div_d * f_c * dt)
                D[t, a, g] = D_val

                dividing_d_cells = prev_D * (1.0 - total_rate_d * dt) * (k_div_d * f_c)
                k1d = k7 + k8 * unrepaired_d
                r_rs_d_val = rs_ub / (1.0 + np.exp(-k1d * g + k0))
                DN[t, 0] += dividing_d_cells * 2.0 * r_rs_d_val
                total_d_sen_rs += dividing_d_cells * 2.0 * r_rs_d_val
                if g + 1 < gen_steps:
                    D[t, 0, g+1] += dividing_d_cells * 2.0 * (1.0 - r_rs_d_val)

        # store bookkeeping
        total_d_sen_ms_array[t] = total_d_sen_ms
        total_d_sen_rs_array[t] = total_d_sen_rs
        total_d_sen_sasp_array[t] = total_d_sen_sasp

    return S, D, SN, DN, SQ, total_d_sen_ms_array, total_d_sen_rs_array, total_d_sen_sasp_array


import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Parameter names and baseline values
# -------------------------

params = {
    "k_div_s": 0.12,
    "k_div_d": 0.1,
    "f_e": 0.0003,
    "p_dif": 0.2,
    "k7": 0.1,
    "k8": 0.05,
    "k0": 5.0,
    "k_a_s": 0.005,
    "k_a_d": 0.01,
    "k_a_sn": 0.01,
    "k_a_dn": 0.05,
    "r_rep_s": 0.5,
    "r_rep_d": 0.4,
    "c_ros": 5.0,
    "k2": 0.00004,
    "k_c": 1.0,
    "resistance_to_apoptosis": 1.0/3.0,
    "r_immuno_surv": 0.9,
    "rs_ub": 0.5,
    "k_q": 0.5,
    "r_ra": 0.001,
    "k_6": 0.02,
    "k_sasp": 0.01,
    "stem_cell_init": 0.02
}




# Explicit bounds for each parameter (lower, baseline, upper)
param_bounds = {
    "k_div_s": [0.012, 0.99],
    "k_div_d": [0.01, 0.99],
    "p_dif": [0.02, 0.99],
    "c_ros": [0.5, 50],
    "f_e": [0.00003, 0.003],
    "k2": [0.000004, 0.0004],
    "k_6": [0.002, 0.2],
    "k7": [0.01, 1],
    "k8": [0.005, 0.5],
    "k0": [0.5, 50],
    "k_a_s": [0.0005, 0.05],
    "k_a_d": [0.001, 0.1],
    "k_a_sn": [0.001, 0.1],
    "k_a_dn": [0.005, 0.5],
    "r_rep_s": [0.05, 0.99],
    "r_rep_d": [0.04, 0.99],
    "k_c": [0.1, 1],
    "resistance_to_apoptosis": [0.033, 0.99],
    "r_immuno_surv": [0.09, 0.99],
    "rs_ub": [0.05, 0.99],
    "k_q": [0.05, 0.99],
    "r_ra": [0.0001, 0.01],
    "k_sasp": [0.001, 0.1],
    "stem_cell_init": [0.002, 0.2]
}



# timestep for 80 years
t_80 = 4000
midline_sen_ratio = None

# store results
sen_ratios = {}
sen_curves = {}

# -------------------------
# Run Baseline
# -------------------------
if __name__ == "__main__":
    t0 = time.time()
    S, D, SN, DN, SQ, total_d_sen_ms, total_d_sen_rs, total_d_sen_sasp = run_simulation_numba(
            initial, int(t_max/dt)+1, age_steps, 80+1, dt, **params)
    t1 = time.time()
    print("Simulation finished in (s):", t1 - t0)
    print("Totals at final timestep:", np.sum(S[4000]), np.sum(D[4000]), np.sum(SQ[4000]), np.sum(SN[4000]), np.sum(DN[4000]))
    sen_total_time = np.sum(SN + DN, axis=1)         
    proliferative_total_time = np.sum(S + D + SQ, axis=(1, 2))  
    total_time = proliferative_total_time + sen_total_time                
    total_time_safe = np.where(total_time == 0, np.nan, total_time)
    sen_curve = sen_total_time / total_time_safe   

    sen_ratios['baseline'] = sen_curve[t_80]
    sen_curves['baseline'] = sen_curve


# -------------------------
# Local sensitivity
# -------------------------
for name, bounds in param_bounds.items():
    low_val, high_val = bounds
    for label, val in zip(["low", "high"], [low_val, high_val]):
        if name == 'k_c' and label =='high':
            continue
        if name =='k7' and label=='high':
            continue
        # copy baseline params and update
        param_copy = params.copy()
        param_copy[name] = val

        # run simulation
        S_mod, D_mod, SN_mod, DN_mod, SQ_mod, *_ = run_simulation_numba(
            initial, int(t_max/dt)+1, age_steps, 80+1, dt,
            **param_copy
        )

        # numerator: total senescent cells at each time (sum over age)
        sen_total_time = np.sum(SN_mod + DN_mod, axis=1)          # shape (time,)
        # denominator: total cells (sum proliferative pools over age+gen, then add senescent totals)
        proliferative_total_time = np.sum(S_mod + D_mod + SQ_mod, axis=(1, 2))  # shape (time,)
        total_time = proliferative_total_time + sen_total_time                 # shape (time,)
        # guard against zeros to avoid divide-by-zero
        total_time_safe = np.where(total_time == 0, np.nan, total_time)
        # senescence ratio curve (time series)
        sen_curve = sen_total_time / total_time_safe   # shape (time,)
        sen_curves[name + label] = sen_curve
        sen_ratio = sen_curve[t_80]
        sen_ratios[name + label] = sen_ratio

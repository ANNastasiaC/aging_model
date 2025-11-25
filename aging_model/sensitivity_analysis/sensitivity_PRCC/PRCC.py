import numpy as np
import pandas as pd
from scipy.stats import qmc
from sklearn.linear_model import LinearRegression

# -----------------------
# 1. LHS sampling
# -----------------------
N = 500  # number of simulations
param_names = list(param_bounds.keys())
d = len(param_names)

sampler = qmc.LatinHypercube(d, seed=42)
lhs_samples = sampler.random(N)

# scale to actual param ranges
param_sets = []
for row in lhs_samples:
    param_dict = {name: low + row[i]*(high-low) for i, (name, (low, high)) in enumerate(param_bounds.items())}
    param_sets.append(param_dict)

# -----------------------
# 2. Run simulations
# -----------------------
Y = []
bad_indices = []

for i, p in enumerate(param_sets):
    try:
        S, D, SN, DN, SQ, *_ = run_simulation_numba(
            initial, time_steps, age_steps, gen_steps, dt, **p
        )
        total_cells = np.sum(S + D + SQ, axis=2) + SN + DN
        y = np.sum(SN[t_80] + DN[t_80]) / np.sum(total_cells[t_80])
        if np.isnan(y) or np.isinf(y):
            raise ValueError("NaN/Inf output")
        Y.append(y)
    except Exception as e:
        bad_indices.append(i)
        Y.append(np.nan)  # placeholder

# -----------------------
# 3. Filter out NaNs
# -----------------------
Y = np.array(Y)
valid_idx = ~np.isnan(Y)
Y = Y[valid_idx]
X = np.array([[ps[name] for name in param_names] for ps in param_sets])
X = X[valid_idx, :]

print(f"Removed {len(bad_indices)} unstable runs:")
for idx in bad_indices:
    print(f"Run {idx}: {param_sets[idx]}")

# -----------------------
# 4. Rank-transform
# -----------------------
X_rank = np.apply_along_axis(lambda col: pd.Series(col).rank().values, 0, X)
Y_rank = pd.Series(Y).rank().values

# -----------------------
# 5. Compute PRCC
# -----------------------
PRCC = {}
for j, name in enumerate(param_names):
    idx_other = [i for i in range(d) if i != j]
    
    model_x = LinearRegression().fit(X_rank[:, idx_other], X_rank[:, j])
    res_x = X_rank[:, j] - model_x.predict(X_rank[:, idx_other])

    model_y = LinearRegression().fit(X_rank[:, idx_other], Y_rank)
    res_y = Y_rank - model_y.predict(X_rank[:, idx_other])

    PRCC[name] = np.corrcoef(res_x, res_y)[0, 1]

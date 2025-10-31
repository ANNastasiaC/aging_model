import math
import numpy as np

initial=1e4

# Time parameters
t_max = 5000  # Timeframe of entire model is 100 years
dt = 1.0     # Time step, 0.02 years
time_steps = int(t_max / dt) + 1  # Include t=0
 
# Age parameters
a_max = 50   # Maximum age of a cell is 1 year
da = 1.0      # Age step, 0.02 years
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation number of cell 
dg = 1        
gen_steps = int(g_max / dg) + 1  

# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)


p_dif = 0.2   # probability of daughter cell being differentiated cell

k_6= 0.02 # probability that a given base pair of damage leads to senescence

#These three paremeters together guarantee that the sigmoid function describing probability of repicative senescence is centered around 50 generations
k7 = 0.1 # k_1 in 2.2.4
k8 = 0.05 # k_2 in 2.2.4
k0 = 5  # k_0 in 2.2.4 

#apoptosis constants (probability that a given base pair of damage would lead to apoptosis in 0.02 years)
k_a_s = 0.005  # for stem
k_a_d = 0.01 # differentiated 
k_a_sn = 0.01 # stem senescent 
k_a_dn = 0.05 # differentiated senescent

# Repair rate adjusted to damage
r_rep_s = 0.5 # for stem cells, 90% of damage are immediately repaired
r_rep_d = 0.4 # for differentiated cells, 99% of damage are immediately repaired

# DNA damage rate for a base pair in a 0.02 years caused by external factors other than ROS
def f_e(t):
    return 0.0003 

c_ros = 5 #ROS concentration
k2 =  0.00004

k_c= 1 #maximum apoptosis rate caused by crowding over 0.02 years

# senescent cells 
resistance_to_apoptosis = 1/3 # senescent cells are resistent to apoptosis caused by DNA damage and crowding, so the apoptosis rate is 1/3 that of other cells
r_immuno_surv=0.9 # immuno surveillance cause 90% of senescent cells to be eliminated in 0.02 years

rs_ub=0.5 # maximum rate of replicative senescence is 0.5

k_q=0.5 # maximum quiescent rate of stem cells in 0.02 years

r_ra=0.001  # reactivation rate of quiescent stem cells in 0.02 years
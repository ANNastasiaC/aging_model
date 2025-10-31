# %%
S_list=[]
D_list=[]
SQ_list=[]
SN_list=[]
DN_list=[]

# %%
# Generate values for t
weeks= 10
d_stem= 0.005



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
# Generate values for t
weeks= 10
d_stem= 0.015



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%


# %%
# Generate values for t
weeks= 10
d_stem= 0.03



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
# Generate values for t
weeks= 50
d_stem= 0.005



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%


# %%


# %%
# Generate values for t
weeks= 50
d_stem= 0.015



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
# Generate values for t
weeks= 50
d_stem= 0.03



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%


# %%
# Generate values for t
weeks= 90
d_stem= 0.005



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
# Generate values for t
weeks= 90
d_stem= 0.015



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%


# %%


# %%
# Generate values for t
weeks= 90
d_stem= 0.03



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
# Generate values for t
weeks= 90
d_stem= 0



import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
initial=1e4
t_values = np.linspace(0, 100, 100)  # From 0 to 15
# Calculate upper limit values using vectorization
def upper_limit(t):
    year=t/50
    weight= 70-80/(1+math.exp((year-7)/5))
    return initial/5*weight

# -*- coding: utf-8 -*- # stem cell depletion not right
"""
Mathematical Modeling of Cellular Aging

This script simulates the dynamics of stem cells, differentiated cells,
and their senescent counterparts over time, age, and generation, based
on the provided mathematical model with corrected parameters.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Define model parameters
# Time parameters
t_max = 5000  # Maximum time
dt = 1.0     # Time step
time_steps = int(t_max / dt) + 1  # Include t=0

# Age parameters
a_max = 50   # Maximum age # max
da = 1.0      # Age step
age_steps = int(a_max / da) + 1  # Include a=0

# Generation parameters
g_max = 80    # Maximum generation
dg = 1        # Generation step
gen_steps = int(g_max / dg) + 1  # Include g=0


# Define grids
t_grid = np.arange(0, t_max + dt, dt)
a_grid = np.arange(0, a_max + da, da)
g_grid = np.arange(0, g_max + dg, dg)

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?


p_dif = 0.2   # probability of daughter cell being differentiated cell?
#r_dm(t)=k2 * c_ros, use k_6 to find r_ms_s(t)

# this is fixed
bp_per_year= 5.1e-10
bp_per_dt= bp_per_year*0.02 # let this be r_dm(t)
c_ros = 5
k2 = bp_per_dt/ c_ros #r damage should be bp_per_dt
total_bp=3e9

k_6= 0.02 # probability that a given bp damage leads to senescence

# so rdm is base pair damage. For each base pair, probability does not lead to r_ms is 1-k_6. prob that no base pair leads to r ms is 

#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

k7 = 0.1 # the telemere thing?
k8 = 0.05 # just scale it to k_7 so that it is reasonable. k7+k8 should be reasonable
k0 = 5  # Offset in the sigmoid function for replicative senescence? #50 generations is 

k_a_s = 0.005 # should it be exponential? what percentage of damage leads to apoptosis?
k_a_d = 0.01
k_a_sn = 0.01
k_a_dn = 0.05


# Repair rate adjusted to damage
r_rep_s = 0.99
r_rep_d = 0.90


# External factors causing DNA damage (adjusted to increase over time)
def f_e(t):
    return 0.0003 #0.0001 works # less in the beginning

# DNA damage rate
def r_dm(t): # in basepair per 0.02
    return max(0,f_e(t) + k2 * c_ros)

# DNA damage-induced senescence rates
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)

k_c=0.1 #apoptosis caused by crowding
def r_a_s(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) # when crowded, f_c is 0. when sparse, f_c is like 1

resistance_to_apoptosis = 1/3 #combined effect of resistance to apoptosis and increased immune surveillance

def r_a_d(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c) )# -*- coding: utf-8 -*-

r_immuno_surv=0.9 # change this part
def r_a_sn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rate parameter
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))


rs_ub=0.5

def r_rs_s(g, t):
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t):
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub

k_q=0.3
def r_q(t,f_c):
    return max(0,k_q*(1-f_c)) # max quescence rate is 0.5
r_ra=0.001
def r_reactivate(t):
    return  r_ra
    # perturbation would be r_ra*(1+math.sin(t))

k_sasp=0.01

def r_s_sasp(sen_ratio):
    return sen_ratio*k_sasp
# Time iteration
accumulated_dna_health=1

for t in range(1, time_steps):
    if 2000<t and 2000+weeks>=t:
        S[t,0,0]+=(np.sum(S[t-1,:,:])+np.sum(SQ[t-1,:,:]))*d_stem
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(accumulated_dna_health)
    total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


S_list.append(S)
D_list.append(D)
SQ_list.append(SQ)
SN_list.append(SN)
DN_list.append(DN)

# %%
label_lst=[
'0.2 years, 0.5% ',
'0.2 years, 1.5% ',
'0.2 years, 3.0% ',
'1.0 years, 0.5%',
'1.0 years, 1.5% ',
'1.0 years, 3.0% ',
'1.8 years, 0.5%',
'1.8 years, 1.5% ',
'1.8 years, 3.0% ',
'control'

]

# %%
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Create sample three-dimensional (A) and two-dimensional (B) arrays
def sen_ration(A,B):

    # Step 2: Sum over the last two dimensions for A
    sum_A = np.sum(A, axis=(1, 2))  # Resulting in a 1D array of shape (4,)

    # Step 3: Sum over the second dimension for B
    sum_B = np.sum(B, axis=1)        # Resulting in a 1D array of shape (4,)

    # Step 4: Calculate the ratio
    # To avoid division by zero, we can use np.where or add a small epsilon
    epsilon = 1e-10
    ratios = sum_B / (sum_A + epsilon)  # Element-wise division
    return ratios

plt.figure(figsize=(6, 4))


for i in range(len(S_list)):
    D_ratios=sen_ration(D_list[i],DN_list[i])
    plt.plot(range(len( D_ratios)),  D_ratios, marker='s', label=label_lst[i], markersize=0, linewidth=0.7)


# Adding titles and labels
#plt.title('Line Plot of Ratios between Summed Arrays')
plt.xlabel('Human age(years)')
plt.ylabel('Senecence ratio')
plt.grid()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

# Show legend
plt.legend(loc='upper left')

# Display the plot
plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt

# Step 1: Create sample three-dimensional (A) and two-dimensional (B) arrays
def sen_ration(A,B):

    # Step 2: Sum over the last two dimensions for A
    sum_A = np.sum(A, axis=(1, 2))  # Resulting in a 1D array of shape (4,)

    # Step 3: Sum over the second dimension for B
    sum_B = np.sum(B, axis=1)        # Resulting in a 1D array of shape (4,)

    # Step 4: Calculate the ratio
    # To avoid division by zero, we can use np.where or add a small epsilon
    epsilon = 1e-10
    ratios = sum_B / (sum_A + epsilon)  # Element-wise division
    return ratios

plt.figure(figsize=(12, 8))


for i in range(len(S_list)-1):
    D_ratios=sen_ration(D_list[i],DN_list[i])
    control_D_ratios=sen_ration(D_list[-1],DN_list[-1])
    plt.plot(range(len( D_ratios)),  D_ratios/control_D_ratios, marker='s', label=label_lst[i], markersize=0, linewidth=0.7)


# Adding titles and labels
#plt.title('Line Plot of Ratios between Summed Arrays')
plt.xlabel('Human age(years)')
plt.ylabel('Relative Senecence Ratio')
plt.grid()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

# Show legend
plt.legend()

# Display the plot
plt.show()



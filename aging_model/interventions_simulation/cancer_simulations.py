# %%
# Generate values for t
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

CS = np.zeros((time_steps, age_steps, gen_steps))
CD = np.zeros((time_steps, age_steps, gen_steps))
CSN = np.zeros((time_steps, age_steps))
CDN = np.zeros((time_steps, age_steps))
CSQ = np.zeros((time_steps, age_steps, gen_steps))


# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells
CS[0, 0, 0] = 5  # Initial number of stem cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?

def r_div_cs(a, f_c):
    return 0.06# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_cd(a, f_c):
    return 0.2 # what ratio of differentiated cells of all age would divide over 0.02 years?


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

    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    CSQ[t,:,:]=CSQ[t-1,:,:]* (1- r_ms_s(current_time - dt)*dt) 
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
        CSN[t,a] = CSN[t-1,a-1] 
        CDN[t,a] = CDN[t-1,a-1] 
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            CS[t, a, g] = CS[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_cs(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*CSQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            CSN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*CS[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*CSQ[t-1,a-1,g]
            CDN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*CD[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            CSQ[t,a,g]+= CS[t - 1, a - 1, g]  * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*CSQ[t-1,a,g]*dt
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_cs_cells= CS[t - 1, a - 1, g] * (1 - total_rate_s * dt) *  r_div_cs(a, f_c)
            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
                CS[t,0,g+1] += dividing_cs_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt))  # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN
            CD[t,0,0] += dividing_cs_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt))

            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            CSN[t,0] += dividing_cs_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            CDN[t,0] += dividing_cs_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
          

            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
            CD[t, a, g] = CD[t - 1, a - 1, g] * (1 - total_rate_d * dt) * (1-r_div_cd(a, f_c)*dt) 
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            dividing_cd_cells= CD[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_cd(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            CDN[t,0] += dividing_cd_cells* 2 *  r_rs_d(g, current_time - dt)
           
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)

            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
                CD[t,0,g+1] += dividing_cd_cells*2*  (1-r_rs_d(g, current_time - dt))
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    #print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(np.sum(CSQ[t,:,:]),np.sum(CS[t,:,:]),np.sum(CD[t,:,:]),np.sum(CSN[t,:]),np.sum(CDN[t,:]))
    #print(accumulated_dna_health)
    #total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


# %%

import numpy as np
import matplotlib.pyplot as plt


# Assuming D is defined elsewhere; for example initialization of D for demonstration
# D = np.random.rand(3001, 81, 100)  # Example initialization of D

# Example 2-layer list representing Z values
Z = np.sum(D[:, :, :], axis=1)

# Create a heatmap using imshow
plt.figure(figsize=(6, 4))  # Adjust the figure size for better readability
plt.imshow(Z, cmap='viridis', aspect='auto', origin='lower')

# Customize the plot

y_ticks = np.arange(0, 5001, 500)
print(y_ticks)
plt.yticks(y_ticks, [int(tick / 50) for tick in y_ticks]) 
plt.colorbar(label='Number of cells')  # Add colorbar to indicate value scale
plt.xlabel('Generation number')
plt.ylabel('Age of human (years)')
#plt.title('Generation number distribution of cells at different times')

# Show the plot
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Example 2-layer list representing Z values
Z = np.sum(D[:,:10,:],axis=2)

# Create grid for X and Y
x = np.arange(Z.shape[1])  # Number of columns
y = np.arange(Z.shape[0])  # Number of rows
X, Y = np.meshgrid(x, y)    # Create a meshgrid

# Create a surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=30, azim=300)
ax.plot_surface(X, Y, Z, cmap='viridis')

# Customize the plot
ax.set_xlabel('Cell age (years)')
ax.set_ylabel('Age of human (years)')
ax.set_zlabel('Number of cells')
#ax.set_title('Surface Plot from 2-Layer List')

y_ticks = ax.get_yticks()  # Get current y-tick locations
ax.set_yticks(y_ticks)  # Set the y-ticks back to the original
ax.set_yticklabels([int(tick / 50) for tick in y_ticks])  # Set new y-tick labels
x_ticks = ax.get_xticks()  # Get current y-tick locations
ax.set_xticks(x_ticks)  # Set the y-ticks back to the original
ax.set_xticklabels([int(tick *0.02) for tick in x_ticks])  # Set new y-tick labels



# Show the plot
plt.show()


# Compute total populations over time
S_total = np.sum(S, axis=(1, 2))
D_total = np.sum(D, axis=(1, 2))
SQ_total = np.sum(SQ, axis=(1, 2))
SN_total = np.sum(SN, axis=1)
DN_total = np.sum(DN, axis=1)
CS_total = np.sum(CS, axis=(1, 2))
CD_total = np.sum(CD, axis=(1, 2))
CSQ_total = np.sum(CSQ, axis=(1, 2))
CSN_total = np.sum(CSN, axis=1)
CDN_total = np.sum(CDN, axis=1)
#print(np.log(DN_total[:200] + 1))
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, S_total, label='Stem Cells (S)')
plt.plot(t_grid, SQ_total, label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, D_total, label='Differentiated Cells (D)')
plt.plot(t_grid, SN_total, label='Senescent Stem Cells (SN)')
plt.plot(t_grid, DN_total, label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, CS_total, label='Cancer Stem Cells (CS)')
plt.plot(t_grid, CSQ_total, label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, CD_total, label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, CSN_total, label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, CDN_total, label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Linear Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
# Change x-ticks to show x/50
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(S_total), label='Stem Cells (S)')
plt.plot(t_grid, np.log(SQ_total), label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, np.log(D_total), label='Differentiated Cells (D)')
plt.plot(t_grid, np.log(SN_total), label='Senescent Stem Cells (SN)')
plt.plot(t_grid, np.log(DN_total), label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, np.log(CS_total), label='Cancer Stem Cells (CS)')
plt.plot(t_grid, np.log(CSQ_total), label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, np.log(CD_total), label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, np.log(CSN_total), label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, np.log(CDN_total), label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Log Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()

senescence_dna_damage = np.zeros(time_steps)
senescence_replicative = np.zeros(time_steps)
senescence_sasp= np.zeros(time_steps)
for t in range(1, time_steps):
    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum
    senescence_dna_damage[t] = np.sum(D[t-1, :,:]) * r_ms_d(t)
    senescence_replicative[t] = 0
    for g in range(1,gen_steps):
         senescence_replicative[t] += np.sum(D[t, :,g]) * 2*r_rs_d(g, t * dt)*r_div_d(a,f_c)*(1-r_ms_d(t-1) +r_s_sasp(sen_ratio) + r_a_d(t-1, f_c))  # Example for each generation
    senescence_sasp[t]= (np.sum(D[t-1, :,:]))*(np.sum(SN[t,:])+np.sum(DN[t,:]))/(np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:]))*k_sasp

plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(senescence_dna_damage), label='Non-telomeric DNA Damage induced senescence', color='purple')
plt.plot(t_grid, np.log(senescence_replicative), label='Replicative senescence', color='brown')
plt.plot(t_grid, np.log(senescence_sasp), label='SASP induced senescence', color='green')
#plt.title('Num')
plt.xlabel('Human age (years)')
plt.ylabel('Number of cells per 0.02 years (log scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.grid()
plt.show()

# %%
# Generate values for t
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

CS = np.zeros((time_steps, age_steps, gen_steps))
CD = np.zeros((time_steps, age_steps, gen_steps))
CSN = np.zeros((time_steps, age_steps))
CDN = np.zeros((time_steps, age_steps))
CSQ = np.zeros((time_steps, age_steps, gen_steps))


# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells
CS[0, 0, 0] = 5  # Initial number of stem cells

def r_div_s(a, f_c):
    return 0.06*f_c# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_d(a, f_c):
    return 0.2*f_c # what ratio of differentiated cells of all age would divide over 0.02 years?

def r_div_cs(a, f_c):
    return 0.06# what ratio of stem cells of all age would divide over 0.02 years? #accounting for quesce 0.5 years (1-prondivide)^25=0.2, divide=1-(0.2)1/25
    
def r_div_cd(a, f_c):
    return 0.2 # what ratio of differentiated cells of all age would divide over 0.02 years?


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

    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t))
    current_time = t * dt

    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum

    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    CSQ[t,:,:]=CSQ[t-1,:,:]
    #print('f_c', f_c)
    #print('r_a_sn(current_time - dt, f_c)',r_a_sn(current_time - dt, f_c))
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
        CSN[t,a] = CSN[t-1,a-1] 
        CDN[t,a] = CDN[t-1,a-1] 
    #print('r_a_d',r_a_d(current_time - dt, f_c))
    total_d_sen_ms=0
    total_d_sen_rs=0
    total_d_sen_sasp=0
    #accumulated_damage+=r_dm(t)/total_bp
    for a in range(1, age_steps):
        for g in range(gen_steps):
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio)
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt
            CS[t, a, g] = CS[t - 1, a - 1, g] *(1-r_div_cs(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*CSQ[t-1,a,g]*dt
            
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 
            #CSN[t,a] += (r_s_sasp(sen_ratio))*CS[t-1,a-1,g]*dt 
            #CDN[t,a] += (r_s_sasp(sen_ratio))*CD[t-1,a-1,g]*dt 
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            CSQ[t,a,g]+= CS[t - 1, a - 1, g]  * (r_q(t, f_c)*dt)- r_reactivate(t)*CSQ[t-1,a,g]*dt
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt


            dividing_cs_cells= CS[t - 1, a - 1, g] *  r_div_cs(a, f_c)
            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) # becomes SN
                CS[t,0,g+1] += dividing_cs_cells* 2 * (1 - p_dif)  # becomes SN
            
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) #becomes DN
            CD[t,0,0] += dividing_cs_cells* 2 *  p_dif 

            #SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            #DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
          

            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            #start of change
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) #delete the r_rs
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
            CD[t, a, g] = CD[t - 1, a - 1, g] * (1-r_div_cd(a, f_c)*dt) 
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)
            dividing_cd_cells= CD[t - 1, a - 1, g] * r_div_cd(a, f_c)
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
           
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)

            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
                CD[t,0,g+1] += dividing_cd_cells*2
            #end of change  
            
            

    print(t)
    #print(r_a_d(t,f_c))
    #print('ms rs sasp',total_d_sen_ms,total_d_sen_rs,total_d_sen_sasp)
    print(np.sum(SQ[t,:,:]),np.sum(S[t,:,:]),np.sum(D[t,:,:]),np.sum(SN[t,:]),np.sum(DN[t,:]))
    print(np.sum(CSQ[t,:,:]),np.sum(CS[t,:,:]),np.sum(CD[t,:,:]),np.sum(CSN[t,:]),np.sum(CDN[t,:]))
    #print(accumulated_dna_health)
    #total=np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:])


# %%


# %%

import numpy as np
import matplotlib.pyplot as plt


# Assuming D is defined elsewhere; for example initialization of D for demonstration
# D = np.random.rand(3001, 81, 100)  # Example initialization of D

# Example 2-layer list representing Z values
Z = np.sum(D[:, :, :], axis=1)

# Create a heatmap using imshow
plt.figure(figsize=(6, 4))  # Adjust the figure size for better readability
plt.imshow(Z, cmap='viridis', aspect='auto', origin='lower')

# Customize the plot

y_ticks = np.arange(0, 5001, 500)
print(y_ticks)
plt.yticks(y_ticks, [int(tick / 50) for tick in y_ticks]) 
plt.colorbar(label='Number of cells')  # Add colorbar to indicate value scale
plt.xlabel('Generation number')
plt.ylabel('Age of human (years)')
#plt.title('Generation number distribution of cells at different times')

# Show the plot
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Example 2-layer list representing Z values
Z = np.sum(D[:,:10,:],axis=2)

# Create grid for X and Y
x = np.arange(Z.shape[1])  # Number of columns
y = np.arange(Z.shape[0])  # Number of rows
X, Y = np.meshgrid(x, y)    # Create a meshgrid

# Create a surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=30, azim=300)
ax.plot_surface(X, Y, Z, cmap='viridis')

# Customize the plot
ax.set_xlabel('Cell age (years)')
ax.set_ylabel('Age of human (years)')
ax.set_zlabel('Number of cells')
#ax.set_title('Surface Plot from 2-Layer List')

y_ticks = ax.get_yticks()  # Get current y-tick locations
ax.set_yticks(y_ticks)  # Set the y-ticks back to the original
ax.set_yticklabels([int(tick / 50) for tick in y_ticks])  # Set new y-tick labels
x_ticks = ax.get_xticks()  # Get current y-tick locations
ax.set_xticks(x_ticks)  # Set the y-ticks back to the original
ax.set_xticklabels([int(tick *0.02) for tick in x_ticks])  # Set new y-tick labels



# Show the plot
plt.show()


# Compute total populations over time
S_total = np.sum(S, axis=(1, 2))
D_total = np.sum(D, axis=(1, 2))
SQ_total = np.sum(SQ, axis=(1, 2))
SN_total = np.sum(SN, axis=1)
DN_total = np.sum(DN, axis=1)
CS_total = np.sum(CS, axis=(1, 2))
CD_total = np.sum(CD, axis=(1, 2))
CSQ_total = np.sum(CSQ, axis=(1, 2))
CSN_total = np.sum(CSN, axis=1)
CDN_total = np.sum(CDN, axis=1)
#print(np.log(DN_total[:200] + 1))
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, S_total, label='Stem Cells (S)')
plt.plot(t_grid, SQ_total, label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, D_total, label='Differentiated Cells (D)')
plt.plot(t_grid, SN_total, label='Senescent Stem Cells (SN)')
plt.plot(t_grid, DN_total, label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, CS_total, label='Cancer Stem Cells (CS)')
plt.plot(t_grid, CSQ_total, label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, CD_total, label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, CSN_total, label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, CDN_total, label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Linear Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
# Change x-ticks to show x/50
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(S_total), label='Stem Cells (S)')
plt.plot(t_grid, np.log(SQ_total), label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, np.log(D_total), label='Differentiated Cells (D)')
plt.plot(t_grid, np.log(SN_total), label='Senescent Stem Cells (SN)')
plt.plot(t_grid, np.log(DN_total), label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, np.log(CS_total), label='Cancer Stem Cells (CS)')
plt.plot(t_grid, np.log(CSQ_total), label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, np.log(CD_total), label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, np.log(CSN_total), label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, np.log(CDN_total), label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Log Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()

senescence_dna_damage = np.zeros(time_steps)
senescence_replicative = np.zeros(time_steps)
senescence_sasp= np.zeros(time_steps)
for t in range(1, time_steps):
    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum
    senescence_dna_damage[t] = np.sum(D[t-1, :,:]) * r_ms_d(t)
    senescence_replicative[t] = 0
    for g in range(1,gen_steps):
         senescence_replicative[t] += np.sum(D[t, :,g]) * 2*r_rs_d(g, t * dt)*r_div_d(a,f_c)*(1-r_ms_d(t-1) +r_s_sasp(sen_ratio) + r_a_d(t-1, f_c))  # Example for each generation
    senescence_sasp[t]= (np.sum(D[t-1, :,:]))*(np.sum(SN[t,:])+np.sum(DN[t,:]))/(np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:]))*k_sasp

plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(senescence_dna_damage), label='Non-telomeric DNA Damage induced senescence', color='purple')
plt.plot(t_grid, np.log(senescence_replicative), label='Replicative senescence', color='brown')
plt.plot(t_grid, np.log(senescence_sasp), label='SASP induced senescence', color='green')
#plt.title('Num')
plt.xlabel('Human age (years)')
plt.ylabel('Number of cells per 0.02 years (log scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.grid()
plt.show()

# %%

nosen_S=S
nosen_D=D
nosen_SN=SN
nosen_DN=DN
nosen_SQ=SQ
import numpy as np
import matplotlib.pyplot as plt

# Assuming D is defined elsewhere; for example initialization of D for demonstration
# D = np.random.rand(3001, 81, 100)  # Example initialization of D

# Example 2-layer list representing Z values
Z = np.sum(D[:, :, :], axis=1)

# Create a heatmap using imshow
plt.figure(figsize=(6, 4))  # Adjust the figure size for better readability
plt.imshow(Z, cmap='viridis', aspect='auto', origin='lower')

# Customize the plot

y_ticks = np.arange(0, 5001, 500)
print(y_ticks)
plt.yticks(y_ticks, [int(tick / 50) for tick in y_ticks]) 
plt.colorbar(label='Number of cells')  # Add colorbar to indicate value scale
plt.xlabel('Generation number')
plt.ylabel('Age of human (years)')
#plt.title('Generation number distribution of cells at different times')

# Show the plot
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Example 2-layer list representing Z values
Z = np.sum(D[:,:10,:],axis=2)

# Create grid for X and Y
x = np.arange(Z.shape[1])  # Number of columns
y = np.arange(Z.shape[0])  # Number of rows
X, Y = np.meshgrid(x, y)    # Create a meshgrid

# Create a surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=30, azim=300)
ax.plot_surface(X, Y, Z, cmap='viridis')

# Customize the plot
ax.set_xlabel('Cell age (years)')
ax.set_ylabel('Age of human (years)')
ax.set_zlabel('Number of cells')
#ax.set_title('Surface Plot from 2-Layer List')

y_ticks = ax.get_yticks()  # Get current y-tick locations
ax.set_yticks(y_ticks)  # Set the y-ticks back to the original
ax.set_yticklabels([int(tick / 50) for tick in y_ticks])  # Set new y-tick labels
x_ticks = ax.get_xticks()  # Get current y-tick locations
ax.set_xticks(x_ticks)  # Set the y-ticks back to the original
ax.set_xticklabels([int(tick *0.02) for tick in x_ticks])  # Set new y-tick labels



# Show the plot
plt.show()


# Compute total populations over time
S_total = np.sum(S, axis=(1, 2))
D_total = np.sum(D, axis=(1, 2))
SQ_total = np.sum(SQ, axis=(1, 2))
SN_total = np.sum(SN, axis=1)
DN_total = np.sum(DN, axis=1)
CS_total = np.sum(CS, axis=(1, 2))
CD_total = np.sum(CD, axis=(1, 2))
CSQ_total = np.sum(CSQ, axis=(1, 2))
CSN_total = np.sum(CSN, axis=1)
CDN_total = np.sum(CDN, axis=1)
#print(np.log(DN_total[:200] + 1))
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, S_total, label='Stem Cells (S)')
plt.plot(t_grid, SQ_total, label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, D_total, label='Differentiated Cells (D)')
plt.plot(t_grid, SN_total, label='Senescent Stem Cells (SN)')
plt.plot(t_grid, DN_total, label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, CS_total, label='Cancer Stem Cells (CS)')
plt.plot(t_grid, CSQ_total, label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, CD_total, label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, CSN_total, label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, CDN_total, label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Linear Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
# Change x-ticks to show x/50
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()
# Plot total populations over time
plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(S_total), label='Stem Cells (S)')
plt.plot(t_grid, np.log(SQ_total), label='Stem Quiescent Cells (SQ)')
plt.plot(t_grid, np.log(D_total), label='Differentiated Cells (D)')
plt.plot(t_grid, np.log(SN_total), label='Senescent Stem Cells (SN)')
plt.plot(t_grid, np.log(DN_total), label='Senescent Differentiated Cells (DN)')
plt.plot(t_grid, np.log(CS_total), label='Cancer Stem Cells (CS)')
plt.plot(t_grid, np.log(CSQ_total), label='Cancer Stem Quiescent Cells (CSQ)')
plt.plot(t_grid, np.log(CD_total), label='Cancer Differentiated Cells (CD)')
plt.plot(t_grid, np.log(CSN_total), label='Cancer Senescent Stem Cells (CSN)')
plt.plot(t_grid, np.log(CDN_total), label='Cancer Senescent Differentiated Cells (CDN)')
plt.xlabel('Human age (years)')
plt.ylabel('Cell Population')
plt.title('Cell Populations Over Time (Log Scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

plt.grid()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.show()

senescence_dna_damage = np.zeros(time_steps)
senescence_replicative = np.zeros(time_steps)
senescence_sasp= np.zeros(time_steps)
for t in range(1, time_steps):
    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio=previous_sen/previous_sum
    f_c=upper_limit(t)/previous_sum
    senescence_dna_damage[t] = np.sum(D[t-1, :,:]) * r_ms_d(t)
    senescence_replicative[t] = 0
    for g in range(1,gen_steps):
         senescence_replicative[t] += np.sum(D[t, :,g]) * 2*r_rs_d(g, t * dt)*r_div_d(a,f_c)*(1-r_ms_d(t-1) +r_s_sasp(sen_ratio) + r_a_d(t-1, f_c))  # Example for each generation
    senescence_sasp[t]= (np.sum(D[t-1, :,:]))*(np.sum(SN[t,:])+np.sum(DN[t,:]))/(np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:]))*k_sasp

plt.figure(figsize=(6, 4))
plt.plot(t_grid, np.log(senescence_dna_damage), label='Non-telomeric DNA Damage induced senescence', color='purple')
plt.plot(t_grid, np.log(senescence_replicative), label='Replicative senescence', color='brown')
plt.plot(t_grid, np.log(senescence_sasp), label='SASP induced senescence', color='green')
#plt.title('Num')
plt.xlabel('Human age (years)')
plt.ylabel('Number of cells per 0.02 years (log scale)')
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), title='Legend')

x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels

plt.grid()
plt.show()



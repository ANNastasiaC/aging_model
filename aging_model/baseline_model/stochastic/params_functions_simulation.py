
import random

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

initial=1e4 #total number of cells of all states at t=0

def upper_limit(t): 
    # returns the baseline total number of cells (maximum number of cells that our initial population could divide into without causing overcrowding), termed N_b(t) in 2.2.1
    year=t/50 
    weight= 70-80/(1+math.exp((year-7)/5)) 
    return initial/5*weight


# Define model parameters

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

# Initialize arrays
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells, assumed to be 2% of total population
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells



def r_div_s(a, f_c):
    return 0.12*f_c*np.random.normal(1, 0.4) # stem cells of all age would divide over 0.02 years 
    
def r_div_d(a, f_c):
    return 0.1*f_c**np.random.normal(1, 0.4) # differentiated cells of all age would divide over 0.02 years


p_dif = 0.2   # probability of daughter cell being differentiated cell


k_6= 0.02 # probability that a given base pair of damage leads to senescence


#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.

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
r_rep_s = 0.5 # for stem cells, % of damage are immediately repaired
r_rep_d = 0.4 # for differentiated cells, % of damage are immediately repaired


# DNA damage rate for a base pair in a 0.02 years caused by external factors other than ROS
def f_e(t):
    return 0.0003 * np.random.normal(1, 0.4)

c_ros = 5 #ROS concentration
k2 =  0.00004

# Total expected DNA damage in base pairs in a 0.02 years
def r_dm(t):   
    return max(0,f_e(t) + k2 * c_ros) #return 0.0005 

# DNA damage-induced senescence rates

# for stem cells
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage) * np.random.normal(1, 0.4)

# for differentiated cells
def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage) *np.random.normal(1, 0.4)

k_c= 1 #maximum apoptosis rate caused by crowding over 0.02 years #0.3 is good
kdamp= 1
# various apoptosis rates
# crowding factor f_c is defined later in the loop. when infinitely crowded, f_c is 0. when infinitely sparse, f_c is 1

def r_a_s(t,f_c): #stem
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c)**kdamp ) *np.random.normal(1, 0.4)####

def r_a_d(t,f_c): #differntiated 
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c)**kdamp ) *np.random.normal(1, 0.4)####

# senescent cells 
resistance_to_apoptosis = 1/3 # senescent cells are resistent to apoptosis caused by DNA damage and crowding, so the apoptosis rate is 1/3 that of other cells
r_immuno_surv=0.9 # immuno surveillance cause 90% of senescent cells to be eliminated in 0.02 years

def r_a_sn(t,f_c): #stem senescent
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c)**kdamp )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health) *np.random.normal(1, 0.4)

def r_a_dn(t,f_c): #differentiated senescent
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage) 
    return max(0, (a_by_damage + k_c*(1-f_c)**kdamp )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health) *np.random.normal(1, 0.4)




rs_ub=0.5 # maximum rate of replicative senescence is 0.5

# Replicative senescence rates calculations, see 2.2.4
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))

def r_rs_s(g, t): # rate of replicative senescence for stem cells
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub *np.random.normal(1, 0.4)

def r_rs_d(g, t): # rate of replicative senescence for differentiated cells
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub *np.random.normal(1, 0.4)

k_q=0.5 # maximum quiescent rate of stem cells in 0.02 years
def r_q(t,f_c): #quiescent rate of stem cells in 0.02 years effected by crowding
    return max(0,k_q*(1-f_c))  *np.random.normal(1, 0.4)

r_ra=0.001  # reactivation rate of quiescent stem cells in 0.02 years
def r_reactivate(t):
    return  r_ra *np.random.normal(1, 0.4)

k_sasp=0.01 # average number of cells experiencing SASP induced senescence in 0.02 years triggered by the SASP of one senescent cells

def r_s_sasp(sen_ratio):#rate of SASP induced senescnece
    return sen_ratio*k_sasp *np.random.normal(1, 0.4)

# Time iteration
accumulated_dna_health=1 # this approxmiates overall dna health that could be used to estimate immune function during immuno surveillance



for t in range(1, time_steps):
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t)) #updates accumulated dna health given current rate of DNA damage before repair.
    
    current_time = t * dt
    #calculating crowding ratio f_c
    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])
    sen_ratio= previous_sen/previous_sum 

    f_c = upper_limit(t)/previous_sum

    

    #updating num of quiescent stem cells with some cells undergoing DNA damage induced senescence
    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 
    
   

    #updating the age of senescent cells with some cells undergoing apoptosis
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt
    if t<1000 and t%50==0:
        print('t=',t)
        print("accumulated DNA health",accumulated_dna_health)
        print('f_c', f_c)
        print("damage induced senescence rate of stem", r_ms_s(current_time - dt))
        print("damage induced senescence rate of differentiated", r_ms_d(current_time - dt))
        print("SASP senescence rate", r_s_sasp(sen_ratio))
        print("apoptosisi stem", r_a_s(current_time - dt, f_c))
        print("apoptosisi differentiated", r_a_d(current_time - dt, f_c))
        print("apoptosis stem senescent", r_a_sn(current_time - dt, f_c))
        print("apoptosis differentiated senescent", r_a_dn(current_time - dt, f_c))
        print("rate of stem cell division", r_div_s(a, f_c))
        print("rate of senescent cell division", r_div_d(a, f_c))
        print("rate of stem cell quiescence",   r_q(t, f_c))
        print("rate of reactivation", r_reactivate(t))
 
    

    #count the total number of
    total_d_sen_ms=0 #differentiated cells undergoing dna-damage induced senescence
    total_d_sen_rs=0 #differentiated cells undergoing replicative senescence
    total_d_sen_sasp=0 #differentiated cells undergoing SASP induced senescence
    #These three quantities are tracked because it is used to plot figure 7
   
    for a in range(1, age_steps):
        for g in range(gen_steps):

            '''
            1 updates the number of active stem cells
            4 updates the number of quiescent stem cells
            8 updates the number of differetiated cells
            2 adds the SASP induced and DNA damage induced senescent cells to the senescent cell count.
            5 calculates the number of dividing stem cells. they divide into 4 kinds of cells. 6 deals with the two kinds of active cells. 7 deals with the two kinds of senescent cells.
            9 calculates the number of dividing differentiated cells. they divide into 2 kinds of cells. 10 deals with active differentiated cells. 11 deal with senescent differentiated cells.

            '''

    

            # 1 updating the age of the stem cells. Some cells are lost because they are undergoing senescence or apoptosis. 
            # senescence here are either DNA damage induced or SASP induced, replicative senescence is accounted during division
            # of the cells not undergoing senescence, some undergo division r_div_s, quiescence r_q. 
            # There is also an inflow of cells by reactivation of quiescent stem cells r_reactivate
            total_rate_s =  r_ms_s(current_time - dt) + r_a_s(current_time - dt, f_c) + r_s_sasp(sen_ratio) #total senescence rate of stem cells
            S[t, a, g] = S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (1-r_div_s(a, f_c)*dt-r_q(t, f_c)*dt) + r_reactivate(t)*SQ[t-1,a,g]*dt

            # stem cells are fluctuating weirdly. So either div rate, quiescence rate or reactivation rate is weird
            
            # 2 adding cells undergoing senescence to the senescent cell count. for stem cells, the quiescent cells undergoing DNA-damage induced senescence are also counted. 
            SN[t,a] += (r_ms_s(current_time - dt)+r_s_sasp(sen_ratio))*S[t-1,a-1,g]*dt + r_ms_s(current_time - dt)*SQ[t-1,a-1,g]
            DN[t,a] += (r_ms_d(current_time - dt)+r_s_sasp(sen_ratio))*D[t-1,a-1,g]*dt 

            # 3 tracking the causes for fig 7
            total_d_sen_ms+=r_ms_d(current_time - dt)*D[t-1,a-1,g]
            total_d_sen_sasp +=r_s_sasp(sen_ratio)*D[t-1,a-1,g]*dt 
           
            # 4 update the number of quiescent cells
            SQ[t,a,g]+= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * (r_q(t, f_c)*dt)- r_reactivate(t)*SQ[t-1,a,g]*dt

            # 5 update the number of dividing stem cells.
            dividing_s_cells= S[t - 1, a - 1, g] * (1 - total_rate_s * dt) * r_div_s(a, f_c)

            # 6 add the newly divided cells back into the count of total cells
            if g+1<81:
                S[t,0,g+1] += dividing_s_cells* 2 * (1 - p_dif) *(1-  r_rs_s(g, current_time - dt)) 
            D[t,0,0] += dividing_s_cells* 2 *  p_dif *(1-  r_rs_s(g, current_time - dt)) 

            # 7 add the newly senescent cells during division into the count of total senescent cells
            SN[t,0] += dividing_s_cells* 2 * (1 - p_dif) *  r_rs_s(g, current_time - dt)
            DN[t,0] += dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)
            total_d_sen_rs+=dividing_s_cells* 2 *p_dif *  r_rs_s(g, current_time - dt)

            # 8 similar as the first equation
            # updating the age of the stem cells. Some cells are lost because they are undergoing senescence or apoptosis. 
            # of the cells not undergoing senescence, some undergo division r_div_s. 
            total_rate_d = (r_ms_d(current_time - dt) +r_s_sasp(sen_ratio) + r_a_d(current_time - dt, f_c)) 
            D[t, a, g] = D[t - 1, a - 1, g] * (1 - total_rate_d * dt)* (1-r_div_d(a, f_c)*dt) 
           
            # 9 calculate the number of dividing cells
            dividing_d_cells= D[t - 1, a - 1, g] * (1 - total_rate_d * dt) * r_div_d(a, f_c)

            # 10 dividing cells that undergo senescence
            DN[t,0] += dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)
            total_d_sen_rs+=dividing_d_cells* 2 *  r_rs_d(g, current_time - dt)

            #11 dividing cells that do not undergo senescence 
            if g+1<81:
                D[t,0,g+1] += dividing_d_cells*2*  (1-r_rs_d(g, current_time - dt))
    

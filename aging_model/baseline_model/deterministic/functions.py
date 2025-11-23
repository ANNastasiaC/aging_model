import math
import numpy as np

def upper_limit(t): 
    # returns the baseline total number of cells (maximum number of cells that our initial population could divide into without causing overcrowding), termed N_b(t) 
    year=t/50 # the time resolution of the model is 0.02 years
    weight= 70-80/(1+math.exp((year-7)/5)) 
    return initial/5*weight #initial birth weight is approximately 5kg

def r_div_s(a, f_c):
    return 0.12*f_c # % of stem cells of all age would divide over 0.02 years
    
def r_div_d(a, f_c):
    return 0.1*f_c # % of differentiated cells of all age would divide over 0.02 years

# Total expected DNA damage in base pairs in a 0.02 years
def r_dm(t):   
    return max(0,f_e(t) + k2 * c_ros) #return 0.0005 

# DNA damage-induced senescence rates

# for stem cells
def r_ms_s(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    return (1- (1-k_6)** unrepaired_damage)

# for differentiated cells
def r_ms_d(t):
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    return (1- (1-k_6)** unrepaired_damage)


def r_a_s(t,f_c): #stem
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_s )** unrepaired_damage)
    return max(0, a_by_damage + k_c*(1-f_c) ) ####

def r_a_d(t,f_c): #differntiated 
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_d)** unrepaired_damage)
    return max(0,a_by_damage + k_c*(1-f_c)) ####

def r_a_sn(t,f_c): #stem senescent
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_s))
    a_by_damage=(1- (1-k_a_sn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)

def r_a_dn(t,f_c): #differentiated senescent
    unrepaired_damage= max(0, r_dm(t) *(1-r_rep_d))
    a_by_damage=(1- (1-k_a_dn)** unrepaired_damage)
    return max(0, (a_by_damage + k_c*(1-f_c) )*resistance_to_apoptosis+r_immuno_surv*accumulated_dna_health)


# Replicative senescence rates calculations, see 2.2.4
def k1_s(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_s))
def k1_d(t):
    return k7 + k8 * max(0, r_dm(t) *(1-r_rep_d))

def r_rs_s(g, t): # rate of replicative senescence for stem cells
    return 1 / (1 + np.exp(-k1_s(t) * g + k0)) *rs_ub

def r_rs_d(g, t): # rate of replicative senescence for differentiated cells
     return 1 / (1 + np.exp(-k1_d(t) * g + k0))*rs_ub


def r_q(t,f_c): #quiescent rate of stem cells in 0.02 years effected by crowding
    return max(0,k_q*(1-f_c)) 


def r_reactivate(t):
    return  r_ra

k_sasp=0.01 # average number of cells experiencing SASP induced senescence in 0.02 years triggered by the SASP of one senescent cells

def r_s_sasp(sen_ratio):#rate of SASP induced senescnece
    return sen_ratio*k_sasp

# Time iteration
accumulated_dna_health=1 # this approxmiates overall dna health that could be used to estimate immune function during immuno surveillance



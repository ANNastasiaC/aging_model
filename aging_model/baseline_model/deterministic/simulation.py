import math
import numpy as np

# Initialize arrays
#S,stem. D, differentiated. SN, stem senescent. DN, differentiated senescent. SQ, stem quescent.
S = np.zeros((time_steps, age_steps, gen_steps))
D = np.zeros((time_steps, age_steps, gen_steps))
SN = np.zeros((time_steps, age_steps))
DN = np.zeros((time_steps, age_steps))
SQ = np.zeros((time_steps, age_steps, gen_steps))

# Initial conditions
initial=1e4
S[0, 0, 0] = initial*0.02  # Initial number of stem cells, assumed to be 2% of total population
D[0, 0, 0] = initial*0.98    # Initial number of differentiated cells

for t in range(1, time_steps):

    #updates accumulated dna health given current rate of DNA damage before repair.
    accumulated_dna_health=accumulated_dna_health*(1-r_dm(t)) 

    current_time = t * dt
    previous_sum=np.sum(S[t-1,:,:])+np.sum(D[t-1,:,:])+np.sum(SN[t-1,:])+np.sum(DN[t-1,:])+np.sum(SQ[t-1,:,:])
    previous_sen=np.sum(SN[t-1,:])+np.sum(DN[t-1,:])

    #calcualte senescence ratio
    sen_ratio= previous_sen/previous_sum 

    #calculate crowding factor
    f_c = upper_limit(t)/previous_sum 

    #updating num of quiescent stem cells with some cells undergoing DNA damage induced senescence
    SQ[t,:,:]=SQ[t-1,:,:]*(1- r_ms_s(current_time - dt)*dt) 

    #updating the age of senescent cells with some cells undergoing apoptosis
    for a in range(1, age_steps):
        SN[t,a] = SN[t-1,a-1] - SN[t-1,a-1]*r_a_sn(current_time - dt, f_c)*dt
        DN[t,a] = DN[t-1,a-1] - DN[t-1,a-1]*r_a_dn(current_time - dt, f_c)*dt

    #count the total number of, for plotting
    total_d_sen_ms=0 #differentiated cells undergoing dna-damage induced senescence
    total_d_sen_rs=0 #differentiated cells undergoing replicative senescence
    total_d_sen_sasp=0 #differentiated cells undergoing SASP induced senescence
   
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
    
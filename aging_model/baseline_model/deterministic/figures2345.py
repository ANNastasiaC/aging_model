from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt

# Simulated population dynamics of cell types throughout the lifespan
fig, axes = plt.subplots(1, 2, figsize=(12,5), dpi=500)  # independent axes
S_total = np.sum(S, axis=(1, 2))
D_total = np.sum(D, axis=(1, 2))
SQ_total = np.sum(SQ, axis=(1, 2))
SN_total = np.sum(SN, axis=1)
DN_total = np.sum(DN, axis=1)
ax = axes[0] # Linear scale
ax.plot(t_grid, S_total, label='Stem Cells (S)')
ax.plot(t_grid, SQ_total, label='Stem Quiescent Cells (SQ)')
ax.plot(t_grid, D_total, label='Differentiated Cells (D)')
ax.plot(t_grid, SN_total, label='Senescent Stem Cells (SN)')
ax.plot(t_grid, DN_total, label='Senescent Differentiated Cells (DN)')
ax.set_xlabel('Human age (years)')
ax.set_ylabel('Cell Population')
ax.set_title('Linear Scale')
ax.grid(True)
ax.set_xticks(x_ticks)
ax.set_xticklabels([int(tick/50) for tick in x_ticks])
ax = axes[1] # Log scale
ax.plot(t_grid, np.log(S_total), label='Stem Cells (S)')
ax.plot(t_grid, np.log(SQ_total), label='Stem Quiescent Cells (SQ)')
ax.plot(t_grid, np.log(D_total), label='Differentiated Cells (D)')
ax.plot(t_grid, np.log(SN_total), label='Senescent Stem Cells (SN)')
ax.plot(t_grid, np.log(DN_total), label='Senescent Differentiated Cells (DN)')
ax.set_xlabel('Human age (years)')
ax.set_title('Log Scale')
ax.grid(True)
ax.set_xticks(x_ticks)
ax.set_xticklabels([int(tick/50) for tick in x_ticks])
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
plt.tight_layout()
plt.show()


# The distribution of cell generation numbers in cells of all states across time points
Z = np.sum(D[:, :, :], axis=1)
plt.figure(figsize=(6, 4),dpi=500)  # Adjust the figure size for better readability
plt.imshow(Z, cmap='viridis', aspect='auto', origin='lower')
y_ticks = np.arange(0, 5001, 500)
print(y_ticks)
plt.yticks(y_ticks, [int(tick / 50) for tick in y_ticks]) 
plt.colorbar(label='Number of cells')  # Add colorbar to indicate value scale
plt.xlabel('Generation number of cell')
plt.ylabel('Age of human (years)')
plt.show()

# Age distribution of differentiated non-senescent cells over time.
Z = np.sum(D[:, :10, :], axis=2)
x = np.arange(Z.shape[1])  # Number of columns
y = np.arange(Z.shape[0])  # Number of rows
X, Y = np.meshgrid(x, y)    # Create a meshgrid
fig = plt.figure(figsize=(8, 5),dpi=500)
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=30, azim=300)
surf = ax.plot_surface(X, Y, Z, cmap='viridis')
fig.colorbar(surf, ax=ax, shrink=0.8, aspect=15, pad=0.15)
ax.set_xlabel('Cell age (years)')
ax.set_ylabel('Age of human (years)')
ax.set_zlabel('Number of cells', labelpad=10)
ax.set_xlim(0, Z.shape[1]-1)
ax.set_ylim(0, Z.shape[0]-1)
ax.set_zlim(0, Z.max())
y_ticks = ax.get_yticks() 
ax.set_yticks(y_ticks)
ax.set_yticklabels([int(tick / 50) for tick in y_ticks])
x_ticks = ax.get_xticks()
x_ticks_half = x_ticks[::2]
ax.set_xticks(x_ticks_half)
ax.set_xticklabels([tick * 0.02 for tick in x_ticks_half])
plt.show()

# Generation number distribution of differentiated non-senescent cells at 100 years.
A, G = np.meshgrid(a_grid, g_grid)
fig = plt.figure(figsize=(8,5), dpi=500) 
D_distribution = D[-1, :15, :80].T  #
ax2 = fig.add_subplot(121, projection='3d')
ax2.view_init(elev=30, azim=340)
surf2 = ax2.plot_surface(A[:80,:15], G[:80,:15], D_distribution, cmap=cm.inferno)
ax2.set_xlabel('Cell Age (years)')
ax2.set_ylabel('Generation Number')
ax2.set_zlabel('Differentiated Cell Density')
x_ticks = ax2.get_xticks()  
ax2.set_xticks(x_ticks)  
ax2.set_xticklabels([round(tick * 0.02, 2) for tick in x_ticks])  #
cax = fig.add_axes([0.59, 0.15, 0.03, 0.7]) 
fig.colorbar(surf2, cax=cax,shrink=0.5, aspect=5)
plt.show()

#Number of new senescent cells by mechanism.
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
    senescence_sasp[t]= (np.sum(D[t-1, :,:]))*(np.sum(SN[t,:])+np.sum(DN[t,:]))/(np.sum(SQ[t,:,:])+np.sum(S[t,:,:])+np.sum(D[t,:,:])+np.sum(SN[t,:])+np.sum(DN[t,:]))*0.01
plt.figure(figsize=(6,4), dpi=500) 
plt.plot(t_grid, np.log(senescence_dna_damage), label='Non-telomeric DNA Damage induced senescence', color='purple')
plt.plot(t_grid, np.log(senescence_replicative), label='Replicative senescence', color='brown')
plt.plot(t_grid, np.log(senescence_sasp), label='SASP induced senescence', color='green')
#plt.title('Num')
plt.xlabel('Human age (years)')
plt.ylabel('Number of cells per 0.02 years (log scale)')
plt.legend()
x_ticks = np.arange(0, 5001, 500)  # Get current x-tick locations
plt.xticks(x_ticks, [int(tick / 50) for tick in x_ticks])  # Set new x-tick labels
plt.grid()
plt.show()


import scipy
import numpy as np
from pylab import *
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)



SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 22
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize

#extracting exmodel from colormap
cmap = cm.get_cmap('plasma', 11)    # PiYG
Plasma=["" for x in range(11)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Plasma[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Reds', 11)    # PiYG
Reds=["" for x in range(11)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Reds[i]=matplotlib.colors.rgb2hex(rgb)

cmap = cm.get_cmap('Blues', 11)    # PiYG
Blues=["" for x in range(11)]

for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        Blues[i]=matplotlib.colors.rgb2hex(rgb)


#Energy Combined per bond

fig = plt.figure(figsize=(3.2,4.5))
ax =fig.add_axes([0.1, 0.14, 0.8, 0.78],frameon=False,xticklabels=[],yticklabels=[])
ax.tick_params(axis='both',direction = 'in',bottom=False, top=False, left=False, right=False)

ax1 = fig.add_axes([0.21, 0.700, 0.785, 0.295], xticklabels=[])#, ylim=(0, 0.6))
ax2 = fig.add_axes([0.21, 0.400, 0.785, 0.295], xticklabels=[])#, ylim=(0, 0.6))
ax3 = fig.add_axes([0.21, 0.100, 0.785, 0.295])#, xticklabels=[])#, ylim=(0, 0.6))

model = np.loadtxt('./P_atm/ene_contri', comments='#')
trans = np.loadtxt('./P_atm/trans_contri', comments='#')

#Ensemble
ax1.plot(model[:,1],(model[:,6]-trans[:,6]-6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color="purple", label=r'$\mathrm{\langle u_{pair} \rangle}$' )
ax1.plot(model[:,1],(trans[:,6])/2.,linewidth="2.5",color='b', label=r'$\mathrm{\langle u_{trans} \rangle}$' )
ax1.plot(model[:,1],(6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color='g', label=r'$\mathrm{\langle u_{equip} \rangle}$' )
ax1.plot(model[:,1],(model[:,6])/2.,linewidth="2.5",color='k', label=r'$\mathrm{\langle u\rangle}$' )
ax1.set_ylabel(r'$\mathrm{\mathbf{\langle u \rangle}\ [kJ/mol]}$',{'fontsize':MEDIUM_SIZE})
ax1.yaxis.set_label_coords(-0.17, 0.5)
ax1.set_xlabel(r'$\mathrm{\mathbf{T}\ [^o C]}$',{'fontsize':MEDIUM_SIZE})
ax1.set_yticks(np.arange(-20, 12, 10))
ax1.set_ylim([-22, 12])
ax1.set_xticks(np.arange(-100, 100, 50))
ax1.set_xlim([-100, 100])
ax1.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax1.legend( frameon=False, fontsize = "10", ncol =2,labelspacing=0.1, loc = "center", bbox_to_anchor=(0.45, 0.48))

#Without u_{equip}rtion kT term
ax2.plot(model[:,1],(model[:,4]-6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color=Reds[10], label=r'$\mathrm{u_{Cage}}$')
ax2.plot(model[:,1],(model[:,3]-6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color=Reds[7], label=r'$\mathrm{u_{pHB}}$' )
ax2.plot(model[:,1],(model[:,5]-6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color=Reds[4], label=r'$\mathrm{u_{vdW}}$' )
ax2.plot(model[:,1],(model[:,6]-6*(model[:,1]+273.15)*8.3145/1000.)/2.,linewidth="2.5",color='k', label=r'$\mathrm{\langle u \rangle}$' )
ax2.set_ylabel(r'$\mathrm{\mathbf{u}\ [kJ/mol]}$',{'fontsize':MEDIUM_SIZE})
ax2.yaxis.set_label_coords(-0.17, 0.5)
ax2.set_xlabel(r'$\mathrm{\mathbf{T}\ [^o C]}$',{'fontsize':MEDIUM_SIZE})
ax2.set_yticks(np.arange(-20, -5, 5))
ax2.set_ylim([-22, -8])
ax2.set_xticks(np.arange(-100, 100, 50))
ax2.set_xlim([-100, 100])
ax2.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)
ax2.legend( frameon=False, fontsize = "10", ncol =2,labelspacing=0.1, loc = "center", bbox_to_anchor=(0.4, 0.5))


#u_{trans}ational contribution
ax3.plot(trans[:,1],(trans[:,4]*1000)/2.,linewidth="2.5",color=Reds[10], label=r'$\mathrm{u_{Cage}}$')
ax3.plot(trans[:,1],(trans[:,3]*1000)/2.,linewidth="2.5",color=Reds[7], label=r'$\mathrm{u_{pHB}}$' )
ax3.plot(trans[:,1],(trans[:,5]*1000)/2.,linewidth="2.5",color=Reds[4], label=r'$\mathrm{u_{vdW}}$' )
ax3.plot(trans[:,1],(trans[:,6]*1000)/2.,linewidth="2.5",color='k', label=r'$\mathrm{\langle u_{trans} \rangle}$' )
ax3.set_ylabel(r'$\mathrm{\mathbf{u}\ [J/mol]}$',{'fontsize':MEDIUM_SIZE})
ax3.yaxis.set_label_coords(-0.17, 0.5)
ax3.set_xlabel(r'$\mathrm{\mathbf{T}\ [^o C]}$',{'fontsize':MEDIUM_SIZE})
ax3.set_yticks(np.arange(0, 150, 50))
ax3.set_ylim([-10, 150])
ax3.set_xticks(np.arange(-100, 100, 50))
ax3.set_xlim([-100, 100])
ax3.tick_params(axis='both',direction = 'in', which= 'both',bottom=1, top=1, left=1,right=1)

plt.savefig("combined_energy_per_bond_wo_kT.png", dpi =600)
plt.show()
plt.close()
